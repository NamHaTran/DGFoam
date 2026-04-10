/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgVtkLagrangeTools.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "basisFunctions.H"
#include "refCoordTransforms.H"

#include <array>
#include <type_traits>

namespace Foam
{
namespace dgVtkLagrange
{

namespace
{

constexpr unsigned char vtkLagrangeTetraType = 71u;
constexpr unsigned char vtkLagrangeHexType   = 72u;
constexpr unsigned char vtkLagrangeWedgeType = 73u;


inline label vtkTriangleIndex
(
    const std::array<label, 3>& bindex,
    label order
)
{
    label index = 0;
    label maxCoord = order;
    label minCoord = 0;
    const label bmin = min(min(bindex[0], bindex[1]), bindex[2]);

    while (bmin > minCoord)
    {
        index += 3*order;
        maxCoord -= 2;
        ++minCoord;
        order -= 3;
    }

    for (label dim = 0; dim < 3; ++dim)
    {
        if (bindex[(dim + 2)%3] == maxCoord)
        {
            return index;
        }

        ++index;
    }

    for (label dim = 0; dim < 3; ++dim)
    {
        if (bindex[(dim + 1)%3] == minCoord)
        {
            return index + bindex[dim] - (minCoord + 1);
        }

        index += maxCoord - (minCoord + 1);
    }

    return index;
}


inline label vtkTetraIndex
(
    const std::array<label, 4>& bindex,
    label order
)
{
    static constexpr label vertexMaxCoords[4] = {3, 0, 1, 2};
    static constexpr label edgeMinCoords[6][2] =
    {
        {1, 2},
        {2, 3},
        {0, 2},
        {0, 1},
        {1, 3},
        {0, 3}
    };
    static constexpr label edgeCountingCoord[6] = {0, 1, 3, 2, 2, 2};
    static constexpr label faceBCoords[4][3] =
    {
        {0, 2, 3},
        {2, 0, 1},
        {2, 1, 3},
        {1, 0, 3}
    };
    static constexpr label faceMinCoord[4] = {1, 3, 0, 2};

    label index = 0;
    label maxCoord = order;
    label minCoord = 0;
    const label bmin =
        min(min(bindex[0], bindex[1]), min(bindex[2], bindex[3]));

    while (bmin > minCoord)
    {
        index += 2*(order*order + 1);
        maxCoord -= 3;
        ++minCoord;
        order -= 4;
    }

    for (label vertex = 0; vertex < 4; ++vertex)
    {
        if (bindex[vertexMaxCoords[vertex]] == maxCoord)
        {
            return index;
        }

        ++index;
    }

    for (label edge = 0; edge < 6; ++edge)
    {
        if
        (
            bindex[edgeMinCoords[edge][0]] == minCoord
         && bindex[edgeMinCoords[edge][1]] == minCoord
        )
        {
            return index + bindex[edgeCountingCoord[edge]] - (minCoord + 1);
        }

        index += maxCoord - (minCoord + 1);
    }

    for (label face = 0; face < 4; ++face)
    {
        if (bindex[faceMinCoord[face]] == minCoord)
        {
            std::array<label, 3> projectedBIndex;

            for (label i = 0; i < 3; ++i)
            {
                projectedBIndex[i] = bindex[faceBCoords[face][i]] - minCoord;
            }

            return index + vtkTriangleIndex(projectedBIndex, order) - 3*order;
        }

        index += ((order + 1)*(order + 2))/2 - 3*order;
    }

    return index;
}


std::vector<label> vtkTetNodemap(const label nPoints)
{
    label order = -1;

    for (label candidate = 1; candidate <= 32; ++candidate)
    {
        if (((candidate + 1)*(candidate + 2)*(candidate + 3))/6 == nPoints)
        {
            order = candidate;
            break;
        }
    }

    if (order < 1)
    {
        FatalErrorInFunction
            << "Cannot deduce tetrahedral order from " << nPoints
            << " points." << exit(FatalError);
    }

    std::vector<label> nodemap(nPoints, -1);
    label pyfrPointI = 0;

    for (label ir = 0; ir <= order; ++ir)
    {
        const label nQ = order + 1 - ir;

        for (label iq = 0; iq < nQ; ++iq)
        {
            const label nP = nQ - iq;

            for (label ip = 0; ip < nP; ++ip)
            {
                const std::array<label, 4> bindex =
                {
                    ip,
                    iq,
                    ir,
                    order - ip - iq - ir
                };

                const label vtkPointI = vtkTetraIndex(bindex, order);
                nodemap[vtkPointI] = pyfrPointI++;
            }
        }
    }

    for (const label pointI : nodemap)
    {
        if (pointI < 0)
        {
            FatalErrorInFunction
                << "Failed to construct a complete VTK tetra nodemap for order "
                << order << exit(FatalError);
        }
    }

    return nodemap;
}


inline vector canonicalEtaFromXi
(
    const dgCellType type,
    const vector& xi
)
{
    const scalar xi2 = xi.y();
    const scalar xi3 = xi.z();

    switch (type)
    {
        case dgCellType::HEX:
            return xi;

        case dgCellType::PRISM:
        {
            const scalar denom = 1.0 - xi3;

            if (mag(denom) <= SMALL)
            {
                return vector(-1.0, xi2, xi3);
            }

            return xiToEtaPrism(xi);
        }

        case dgCellType::TET:
        {
            const scalar denomEta2 = 1.0 - xi3;
            const scalar denomEta1 = -xi2 - xi3;

            if (mag(denomEta2) <= SMALL)
            {
                return vector(-1.0, -1.0, xi3);
            }

            const scalar eta2 = 2.0*(xi2 + 1.0)/denomEta2 - 1.0;

            // Along the collapsed C-D edge, eta1 is not unique in xi-space,
            // but eta2 and eta3 remain well defined.
            if (mag(denomEta1) <= SMALL)
            {
                return vector(-1.0, eta2, xi3);
            }

            return xiToEtaTet(xi);
        }

        default:
            FatalErrorInFunction
                << "No xi->eta transform is implemented for "
                << cellTypeName(type) << exit(FatalError);
    }

    return vector::zero;
}


inline void writePointComponentTriples
(
    Ostream& os,
    const std::vector<point>& values
)
{
    for (const point& pt : values)
    {
        os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << ' ';
    }

    os << nl;
}


inline void writeVectorComponentTriples
(
    Ostream& os,
    const std::vector<vector>& values
)
{
    for (const vector& v : values)
    {
        os << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
    }

    os << nl;
}


template<class ScalarLike>
void writeScalarValues
(
    Ostream& os,
    const std::vector<ScalarLike>& values
)
{
    for (const ScalarLike& v : values)
    {
        os << v << ' ';
    }

    os << nl;
}


template<class IntType>
void writeIntDataArray
(
    Ostream& os,
    const word& name,
    const word& vtkType,
    const std::vector<IntType>& values
)
{
    os  << "        <DataArray type=\"" << vtkType
        << "\" Name=\"" << name
        << "\" NumberOfComponents=\"1\" format=\"ascii\">" << nl
        << "          ";
    writeScalarValues(os, values);
    os  << "        </DataArray>" << nl;
}


void writeScalarPointDataArray
(
    Ostream& os,
    const ScalarPointData& array
)
{
    os  << "        <DataArray type=\"" << vtkFloatTypeName()
        << "\" Name=\"" << array.name
        << "\" NumberOfComponents=\"1\" format=\"ascii\">" << nl
        << "          ";
    writeScalarValues(os, *array.values);
    os  << "        </DataArray>" << nl;
}


void writeVectorPointDataArray
(
    Ostream& os,
    const VectorPointData& array
)
{
    os  << "        <DataArray type=\"" << vtkFloatTypeName()
        << "\" Name=\"" << array.name
        << "\" NumberOfComponents=\"3\" format=\"ascii\">" << nl
        << "          ";
    writeVectorComponentTriples(os, *array.values);
    os  << "        </DataArray>" << nl;
}

} // End anonymous namespace


const char* vtkFloatTypeName()
{
    return std::is_same<scalar, double>::value ? "Float64" : "Float32";
}


word cellTypeName(const dgCellType type)
{
    switch (type)
    {
        case dgCellType::HEX:
            return "hex";

        case dgCellType::PRISM:
            return "prism";

        case dgCellType::TET:
            return "tet";

        case dgCellType::PYRAMID:
            return "pyramid";

        default:
            return "invalid";
    }
}


bool supportsNativeLagrange(const dgCellType type)
{
    return
        type == dgCellType::HEX
     || type == dgCellType::PRISM
     || type == dgCellType::TET;
}


unsigned char vtkLagrangeType(const dgCellType type)
{
    switch (type)
    {
        case dgCellType::HEX:
            return vtkLagrangeHexType;

        case dgCellType::PRISM:
            return vtkLagrangeWedgeType;

        case dgCellType::TET:
            return vtkLagrangeTetraType;

        default:
            FatalErrorInFunction
                << "No native VTK Lagrange cell code for DG cell type "
                << cellTypeName(type) << exit(FatalError);
    }

    return 0u;
}


label exportOrder(const label pOrder)
{
    return max(label(1), pOrder);
}


std::vector<label> vtkNodemap(const dgCellType type, const label nPoints)
{
    switch (type)
    {
        case dgCellType::HEX:
        {
            switch (nPoints)
            {
                case 8:   return {0, 1, 3, 2, 4, 5, 7, 6};
                case 27:  return {0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15, 12, 14, 10, 16, 4, 22, 13};
                case 64:  return {0, 3, 15, 12, 48, 51, 63, 60, 1, 2, 7, 11, 13, 14, 4, 8, 49, 50, 55, 59, 61, 62, 52, 56, 16, 32, 19, 35, 31, 47, 28, 44, 20, 24, 36, 40, 23, 27, 39, 43, 17, 18, 33, 34, 29, 30, 45, 46, 5, 6, 9, 10, 53, 54, 57, 58, 21, 22, 25, 26, 37, 38, 41, 42};
                case 125: return {0, 4, 24, 20, 100, 104, 124, 120, 1, 2, 3, 9, 14, 19, 21, 22, 23, 5, 10, 15, 101, 102, 103, 109, 114, 119, 121, 122, 123, 105, 110, 115, 25, 50, 75, 29, 54, 79, 49, 74, 99, 45, 70, 95, 30, 35, 40, 55, 60, 65, 80, 85, 90, 34, 39, 44, 59, 64, 69, 84, 89, 94, 26, 27, 28, 51, 52, 53, 76, 77, 78, 46, 47, 48, 71, 72, 73, 96, 97, 98, 6, 7, 8, 11, 12, 13, 16, 17, 18, 106, 107, 108, 111, 112, 113, 116, 117, 118, 31, 32, 33, 36, 37, 38, 41, 42, 43, 56, 57, 58, 61, 62, 63, 66, 67, 68, 81, 82, 83, 86, 87, 88, 91, 92, 93};
                case 216: return {0, 5, 35, 30, 180, 185, 215, 210, 1, 2, 3, 4, 11, 17, 23, 29, 31, 32, 33, 34, 6, 12, 18, 24, 181, 182, 183, 184, 191, 197, 203, 209, 211, 212, 213, 214, 186, 192, 198, 204, 36, 72, 108, 144, 41, 77, 113, 149, 71, 107, 143, 179, 66, 102, 138, 174, 42, 48, 54, 60, 78, 84, 90, 96, 114, 120, 126, 132, 150, 156, 162, 168, 47, 53, 59, 65, 83, 89, 95, 101, 119, 125, 131, 137, 155, 161, 167, 173, 37, 38, 39, 40, 73, 74, 75, 76, 109, 110, 111, 112, 145, 146, 147, 148, 67, 68, 69, 70, 103, 104, 105, 106, 139, 140, 141, 142, 175, 176, 177, 178, 7, 8, 9, 10, 13, 14, 15, 16, 19, 20, 21, 22, 25, 26, 27, 28, 187, 188, 189, 190, 193, 194, 195, 196, 199, 200, 201, 202, 205, 206, 207, 208, 43, 44, 45, 46, 49, 50, 51, 52, 55, 56, 57, 58, 61, 62, 63, 64, 79, 80, 81, 82, 85, 86, 87, 88, 91, 92, 93, 94, 97, 98, 99, 100, 115, 116, 117, 118, 121, 122, 123, 124, 127, 128, 129, 130, 133, 134, 135, 136, 151, 152, 153, 154, 157, 158, 159, 160, 163, 164, 165, 166, 169, 170, 171, 172};
                case 343: return {0, 6, 48, 42, 294, 300, 342, 336, 1, 2, 3, 4, 5, 13, 20, 27, 34, 41, 43, 44, 45, 46, 47, 7, 14, 21, 28, 35, 295, 296, 297, 298, 299, 307, 314, 321, 328, 335, 337, 338, 339, 340, 341, 301, 308, 315, 322, 329, 49, 98, 147, 196, 245, 55, 104, 153, 202, 251, 97, 146, 195, 244, 293, 91, 140, 189, 238, 287, 56, 63, 70, 77, 84, 105, 112, 119, 126, 133, 154, 161, 168, 175, 182, 203, 210, 217, 224, 231, 252, 259, 266, 273, 280, 62, 69, 76, 83, 90, 111, 118, 125, 132, 139, 160, 167, 174, 181, 188, 209, 216, 223, 230, 237, 258, 265, 272, 279, 286, 50, 51, 52, 53, 54, 99, 100, 101, 102, 103, 148, 149, 150, 151, 152, 197, 198, 199, 200, 201, 246, 247, 248, 249, 250, 92, 93, 94, 95, 96, 141, 142, 143, 144, 145, 190, 191, 192, 193, 194, 239, 240, 241, 242, 243, 288, 289, 290, 291, 292, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 29, 30, 31, 32, 33, 36, 37, 38, 39, 40, 302, 303, 304, 305, 306, 309, 310, 311, 312, 313, 316, 317, 318, 319, 320, 323, 324, 325, 326, 327, 330, 331, 332, 333, 334, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 71, 72, 73, 74, 75, 78, 79, 80, 81, 82, 85, 86, 87, 88, 89, 106, 107, 108, 109, 110, 113, 114, 115, 116, 117, 120, 121, 122, 123, 124, 127, 128, 129, 130, 131, 134, 135, 136, 137, 138, 155, 156, 157, 158, 159, 162, 163, 164, 165, 166, 169, 170, 171, 172, 173, 176, 177, 178, 179, 180, 183, 184, 185, 186, 187, 204, 205, 206, 207, 208, 211, 212, 213, 214, 215, 218, 219, 220, 221, 222, 225, 226, 227, 228, 229, 232, 233, 234, 235, 236, 253, 254, 255, 256, 257, 260, 261, 262, 263, 264, 267, 268, 269, 270, 271, 274, 275, 276, 277, 278, 281, 282, 283, 284, 285};
                default:
                    break;
            }

            break;
        }

        case dgCellType::TET:
        {
            return vtkTetNodemap(nPoints);
        }

        case dgCellType::PRISM:
        {
            switch (nPoints)
            {
                case 6:   return {0, 1, 2, 3, 4, 5};
                case 18:  return {0, 2, 5, 12, 14, 17, 1, 4, 3, 13, 16, 15, 6, 8, 11, 7, 10, 9};
                case 40:  return {0, 3, 9, 30, 33, 39, 1, 2, 6, 8, 7, 4, 31, 32, 36, 38, 37, 34, 10, 20, 13, 23, 19, 29, 5, 35, 11, 12, 21, 22, 16, 18, 26, 28, 17, 14, 27, 24, 15, 25};
                case 75:  return {0, 4, 14, 60, 64, 74, 1, 2, 3, 8, 11, 13, 12, 9, 5, 61, 62, 63, 68, 71, 73, 72, 69, 65, 15, 30, 45, 19, 34, 49, 29, 44, 59, 6, 7, 10, 66, 67, 70, 16, 17, 18, 31, 32, 33, 46, 47, 48, 23, 26, 28, 38, 41, 43, 53, 56, 58, 27, 24, 20, 42, 39, 35, 57, 54, 50, 21, 22, 25, 36, 37, 40, 51, 52, 55};
                case 126: return {0, 5, 20, 105, 110, 125, 1, 2, 3, 4, 10, 14, 17, 19, 18, 15, 11, 6, 106, 107, 108, 109, 115, 119, 122, 124, 123, 120, 116, 111, 21, 42, 63, 84, 26, 47, 68, 89, 41, 62, 83, 104, 7, 8, 9, 12, 13, 16, 112, 113, 114, 117, 118, 121, 22, 23, 24, 25, 43, 44, 45, 46, 64, 65, 66, 67, 85, 86, 87, 88, 31, 35, 38, 40, 52, 56, 59, 61, 73, 77, 80, 82, 94, 98, 101, 103, 39, 36, 32, 27, 60, 57, 53, 48, 81, 78, 74, 69, 102, 99, 95, 90, 28, 29, 30, 33, 34, 37, 49, 50, 51, 54, 55, 58, 70, 71, 72, 75, 76, 79, 91, 92, 93, 96, 97, 100};
                case 196: return {0, 6, 27, 168, 174, 195, 1, 2, 3, 4, 5, 12, 17, 21, 24, 26, 25, 22, 18, 13, 7, 169, 170, 171, 172, 173, 180, 185, 189, 192, 194, 193, 190, 186, 181, 175, 28, 56, 84, 112, 140, 34, 62, 90, 118, 146, 55, 83, 111, 139, 167, 8, 9, 10, 11, 14, 15, 16, 19, 20, 23, 176, 177, 178, 179, 182, 183, 184, 187, 188, 191, 29, 30, 31, 32, 33, 57, 58, 59, 60, 61, 85, 86, 87, 88, 89, 113, 114, 115, 116, 117, 141, 142, 143, 144, 145, 40, 45, 49, 52, 54, 68, 73, 77, 80, 82, 96, 101, 105, 108, 110, 124, 129, 133, 136, 138, 152, 157, 161, 164, 166, 53, 50, 46, 41, 35, 81, 78, 74, 69, 63, 109, 106, 102, 97, 91, 137, 134, 130, 125, 119, 165, 162, 158, 153, 147, 36, 37, 38, 39, 42, 43, 44, 47, 48, 51, 64, 65, 66, 67, 70, 71, 72, 75, 76, 79, 92, 93, 94, 95, 98, 99, 100, 103, 104, 107, 120, 121, 122, 123, 126, 127, 128, 131, 132, 135, 148, 149, 150, 151, 154, 155, 156, 159, 160, 163};
                case 288: return {0, 7, 35, 252, 259, 287, 1, 2, 3, 4, 5, 6, 14, 20, 25, 29, 32, 34, 33, 30, 26, 21, 15, 8, 253, 254, 255, 256, 257, 258, 266, 272, 277, 281, 284, 286, 285, 282, 278, 273, 267, 260, 36, 72, 108, 144, 180, 216, 43, 79, 115, 151, 187, 223, 71, 107, 143, 179, 215, 251, 9, 10, 11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 27, 28, 31, 261, 262, 263, 264, 265, 268, 269, 270, 271, 274, 275, 276, 279, 280, 283, 37, 38, 39, 40, 41, 42, 73, 74, 75, 76, 77, 78, 109, 110, 111, 112, 113, 114, 145, 146, 147, 148, 149, 150, 181, 182, 183, 184, 185, 186, 217, 218, 219, 220, 221, 222, 50, 56, 61, 65, 68, 70, 86, 92, 97, 101, 104, 106, 122, 128, 133, 137, 140, 142, 158, 164, 169, 173, 176, 178, 194, 200, 205, 209, 212, 214, 230, 236, 241, 245, 248, 250, 69, 66, 62, 57, 51, 44, 105, 102, 98, 93, 87, 80, 141, 138, 134, 129, 123, 116, 177, 174, 170, 165, 159, 152, 213, 210, 206, 201, 195, 188, 249, 246, 242, 237, 231, 224, 45, 46, 47, 48, 49, 52, 53, 54, 55, 58, 59, 60, 63, 64, 67, 81, 82, 83, 84, 85, 88, 89, 90, 91, 94, 95, 96, 99, 100, 103, 117, 118, 119, 120, 121, 124, 125, 126, 127, 130, 131, 132, 135, 136, 139, 153, 154, 155, 156, 157, 160, 161, 162, 163, 166, 167, 168, 171, 172, 175, 189, 190, 191, 192, 193, 196, 197, 198, 199, 202, 203, 204, 207, 208, 211, 225, 226, 227, 228, 229, 232, 233, 234, 235, 238, 239, 240, 243, 244, 247};
                case 405: return {0, 8, 44, 360, 368, 404, 1, 2, 3, 4, 5, 6, 7, 16, 23, 29, 34, 38, 41, 43, 42, 39, 35, 30, 24, 17, 9, 361, 362, 363, 364, 365, 366, 367, 376, 383, 389, 394, 398, 401, 403, 402, 399, 395, 390, 384, 377, 369, 45, 90, 135, 180, 225, 270, 315, 53, 98, 143, 188, 233, 278, 323, 89, 134, 179, 224, 269, 314, 359, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 25, 26, 27, 28, 31, 32, 33, 36, 37, 40, 370, 371, 372, 373, 374, 375, 378, 379, 380, 381, 382, 385, 386, 387, 388, 391, 392, 393, 396, 397, 400, 46, 47, 48, 49, 50, 51, 52, 91, 92, 93, 94, 95, 96, 97, 136, 137, 138, 139, 140, 141, 142, 181, 182, 183, 184, 185, 186, 187, 226, 227, 228, 229, 230, 231, 232, 271, 272, 273, 274, 275, 276, 277, 316, 317, 318, 319, 320, 321, 322, 61, 68, 74, 79, 83, 86, 88, 106, 113, 119, 124, 128, 131, 133, 151, 158, 164, 169, 173, 176, 178, 196, 203, 209, 214, 218, 221, 223, 241, 248, 254, 259, 263, 266, 268, 286, 293, 299, 304, 308, 311, 313, 331, 338, 344, 349, 353, 356, 358, 87, 84, 80, 75, 69, 62, 54, 132, 129, 125, 120, 114, 107, 99, 177, 174, 170, 165, 159, 152, 144, 222, 219, 215, 210, 204, 197, 189, 267, 264, 260, 255, 249, 242, 234, 312, 309, 305, 300, 294, 287, 279, 357, 354, 350, 345, 339, 332, 324, 55, 56, 57, 58, 59, 60, 63, 64, 65, 66, 67, 70, 71, 72, 73, 76, 77, 78, 81, 82, 85, 100, 101, 102, 103, 104, 105, 108, 109, 110, 111, 112, 115, 116, 117, 118, 121, 122, 123, 126, 127, 130, 145, 146, 147, 148, 149, 150, 153, 154, 155, 156, 157, 160, 161, 162, 163, 166, 167, 168, 171, 172, 175, 190, 191, 192, 193, 194, 195, 198, 199, 200, 201, 202, 205, 206, 207, 208, 211, 212, 213, 216, 217, 220, 235, 236, 237, 238, 239, 240, 243, 244, 245, 246, 247, 250, 251, 252, 253, 256, 257, 258, 261, 262, 265, 280, 281, 282, 283, 284, 285, 288, 289, 290, 291, 292, 295, 296, 297, 298, 301, 302, 303, 306, 307, 310, 325, 326, 327, 328, 329, 330, 333, 334, 335, 336, 337, 340, 341, 342, 343, 346, 347, 348, 351, 352, 355};
                default:
                    break;
            }

            break;
        }

        default:
            break;
    }

    FatalErrorInFunction
        << "No VTK nodemap is available for cell type " << cellTypeName(type)
        << " with " << nPoints << " nodes. "
        << "Current native exporter supports high-order orders 1..8 for "
        << "tet/prism/hex." << exit(FatalError);

    return {};
}


std::vector<vector> pyfrStdElementPoints(const dgCellType type, const label order)
{
    std::vector<vector> pts;

    const label n1 = order + 1;
    const scalar delta = 2.0/max(label(1), order);

    switch (type)
    {
        case dgCellType::HEX:
        {
            pts.reserve(n1*n1*n1);

            for (label k = 0; k < n1; ++k)
            {
                const scalar z = -1.0 + k*delta;

                for (label j = 0; j < n1; ++j)
                {
                    const scalar y = -1.0 + j*delta;

                    for (label i = 0; i < n1; ++i)
                    {
                        const scalar x = -1.0 + i*delta;
                        pts.emplace_back(x, y, z);
                    }
                }
            }

            break;
        }

        case dgCellType::TET:
        {
            const label nPts = (n1*(n1 + 1)*(n1 + 2))/6;
            pts.reserve(nPts);

            for (label ir = 0; ir < n1; ++ir)
            {
                const scalar r = -1.0 + ir*delta;
                const label nQ = n1 - ir;

                for (label iq = 0; iq < nQ; ++iq)
                {
                    const scalar q = -1.0 + iq*delta;
                    const label nP = nQ - iq;

                    for (label ip = 0; ip < nP; ++ip)
                    {
                        const scalar p = -1.0 + ip*delta;
                        pts.emplace_back(p, q, r);
                    }
                }
            }

            break;
        }

        case dgCellType::PRISM:
        {
            const label nTri = n1*(n1 + 1)/2;
            pts.reserve(n1*nTri);

            for (label ir = 0; ir < n1; ++ir)
            {
                const scalar r = -1.0 + ir*delta;

                for (label iq = 0; iq < n1; ++iq)
                {
                    const scalar q = -1.0 + iq*delta;
                    const label nP = n1 - iq;

                    for (label ip = 0; ip < nP; ++ip)
                    {
                        const scalar p = -1.0 + ip*delta;
                        pts.emplace_back(p, q, r);
                    }
                }
            }

            break;
        }

        default:
            FatalErrorInFunction
                << "No PyFR-style standard-element points are implemented for "
                << cellTypeName(type) << exit(FatalError);
    }

    return pts;
}


vector xiToEta(const dgCellType type, const vector& xi)
{
    return canonicalEtaFromXi(type, xi);
}


void computeBasisAt
(
    const dgCellType type,
    const vector& eta,
    const label pOrder,
    List<scalar>& basis
)
{
    List<scalar> dEta1;
    List<scalar> dEta2;
    List<scalar> dEta3;

    switch (type)
    {
        case dgCellType::HEX:
            computeHexBasisAndDerivatives
            (
                eta,
                pOrder,
                basis,
                dEta1,
                dEta2,
                dEta3
            );
            break;

        case dgCellType::PRISM:
            computePrismBasisAndDerivatives
            (
                eta,
                pOrder,
                basis,
                dEta1,
                dEta2,
                dEta3
            );
            break;

        case dgCellType::TET:
            computeTetBasisAndDerivatives
            (
                eta,
                pOrder,
                basis,
                dEta1,
                dEta2,
                dEta3
            );
            break;

        default:
            FatalErrorInFunction
                << "Unsupported basis evaluation for cell type "
                << cellTypeName(type) << exit(FatalError);
    }
}


void writeVtuPiece
(
    const fileName& filePath,
    const scalar timeValue,
    const GridBuffers& grid,
    const std::vector<ScalarPointData>& scalarArrays,
    const std::vector<VectorPointData>& vectorArrays
)
{
    mkDir(filePath.path());

    OFstream os(filePath);
    os.precision(16);

    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"UnstructuredGrid\" version=\"2.1\" "
        << "byte_order=\"LittleEndian\">" << nl
        << "  <UnstructuredGrid>" << nl
        << "    <FieldData>" << nl
        << "      <DataArray type=\"Float64\" Name=\"TimeValue\" "
        << "NumberOfTuples=\"1\" NumberOfComponents=\"1\" format=\"ascii\">"
        << nl
        << "        " << timeValue << nl
        << "      </DataArray>" << nl
        << "    </FieldData>" << nl
        << "    <Piece NumberOfPoints=\"" << grid.points.size()
        << "\" NumberOfCells=\"" << grid.cellTypes.size() << "\">" << nl
        << "      <Points>" << nl
        << "        <DataArray type=\"" << vtkFloatTypeName()
        << "\" NumberOfComponents=\"3\" format=\"ascii\">" << nl
        << "          ";

    writePointComponentTriples(os, grid.points);

    os  << "        </DataArray>" << nl
        << "      </Points>" << nl
        << "      <Cells>" << nl;

    writeIntDataArray(os, "connectivity", "Int64", grid.connectivity);
    writeIntDataArray(os, "offsets", "Int64", grid.offsets);

    os  << "        <DataArray type=\"UInt8\" Name=\"types\" "
        << "NumberOfComponents=\"1\" format=\"ascii\">" << nl
        << "          ";

    for (const unsigned char cellType : grid.cellTypes)
    {
        os << label(cellType) << ' ';
    }

    os << nl;
    os  << "        </DataArray>" << nl
        << "      </Cells>" << nl
        << "      <CellData>" << nl;

    writeIntDataArray(os, "Partition", "Int32", grid.partitions);

    os  << "      </CellData>" << nl
        << "      <PointData>" << nl;

    for (const ScalarPointData& array : scalarArrays)
    {
        writeScalarPointDataArray(os, array);
    }

    for (const VectorPointData& array : vectorArrays)
    {
        writeVectorPointDataArray(os, array);
    }

    os  << "      </PointData>" << nl
        << "    </Piece>" << nl
        << "  </UnstructuredGrid>" << nl
        << "</VTKFile>" << nl;
}


void writePvtu
(
    const fileName& filePath,
    const scalar timeValue,
    const label nProcs,
    const word& pieceStem,
    const std::vector<PointDataSchema>& pointData
)
{
    mkDir(filePath.path());

    OFstream os(filePath);

    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"PUnstructuredGrid\" version=\"2.1\" "
        << "byte_order=\"LittleEndian\">" << nl
        << "  <PUnstructuredGrid GhostLevel=\"0\">" << nl
        << "    <FieldData>" << nl
        << "      <DataArray type=\"Float64\" Name=\"TimeValue\" "
        << "NumberOfTuples=\"1\" NumberOfComponents=\"1\" format=\"ascii\">"
        << nl
        << "        " << timeValue << nl
        << "      </DataArray>" << nl
        << "    </FieldData>" << nl
        << "    <PPoints>" << nl
        << "      <PDataArray type=\"" << vtkFloatTypeName()
        << "\" NumberOfComponents=\"3\"/>" << nl
        << "    </PPoints>" << nl
        << "    <PCells>" << nl
        << "      <PDataArray type=\"Int64\" Name=\"connectivity\"/>" << nl
        << "      <PDataArray type=\"Int64\" Name=\"offsets\"/>" << nl
        << "      <PDataArray type=\"UInt8\" Name=\"types\"/>" << nl
        << "    </PCells>" << nl
        << "    <PCellData>" << nl
        << "      <PDataArray type=\"Int32\" Name=\"Partition\" "
        << "NumberOfComponents=\"1\"/>" << nl
        << "    </PCellData>" << nl
        << "    <PPointData>" << nl;

    for (const PointDataSchema& array : pointData)
    {
        os  << "      <PDataArray type=\"" << vtkFloatTypeName()
            << "\" Name=\"" << array.name
            << "\" NumberOfComponents=\"" << array.nComponents
            << "\"/>" << nl;
    }

    os  << "    </PPointData>" << nl;

    for (label procI = 0; procI < nProcs; ++procI)
    {
        os  << "    <Piece Source=\""
            << pieceStem << "_" << procI << ".vtu\"/>" << nl;
    }

    os  << "  </PUnstructuredGrid>" << nl
        << "</VTKFile>" << nl;
}

} // End namespace dgVtkLagrange
} // End namespace Foam

// ************************************************************************* //
