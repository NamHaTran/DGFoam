/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgEquiSpacedTessellation.H"

#include "DynamicList.H"
#include "Jacobian.H"
#include "basisFunctions.H"
#include "refCoordTransforms.H"

#include <initializer_list>

namespace Foam
{

namespace
{

constexpr unsigned char vtkTetType = 10u;
constexpr unsigned char vtkHexType = 12u;
constexpr unsigned char vtkWedgeType = 13u;
constexpr unsigned char vtkPyramidType = 14u;
constexpr unsigned char vtkTriType = 5u;
constexpr unsigned char vtkQuadType = 9u;

inline label trianglePointCount(const label m)
{
    return (m + 1)*(m + 2)/2;
}

inline label triangleRowOffset(const label m, const label row)
{
    return row*(m + 1) - row*(row - 1)/2;
}

inline label prismPointIndex
(
    const label n,
    const label layer,
    const label i,
    const label j
)
{
    return layer*trianglePointCount(n) + triangleRowOffset(n, j) + i;
}

inline label pyramidPointIndex
(
    const label n,
    const label layer,
    const label i,
    const label j
)
{
    label offset = 0;

    for (label k = 0; k < layer; ++k)
    {
        const label width = n - k + 1;
        offset += width*width;
    }

    const label width = n - layer + 1;
    return offset + j*width + i;
}

inline label tetPointIndex
(
    const label n,
    const label layer,
    const label i,
    const label j
)
{
    label offset = 0;

    for (label k = 0; k < layer; ++k)
    {
        offset += trianglePointCount(n - k);
    }

    const label m = n - layer;
    return offset + triangleRowOffset(m, j) + i;
}

inline void appendSubCell
(
    DynamicList<dgEquiSpacedTessellation::VtkSubCell>& cells,
    std::initializer_list<label> localVerts,
    const unsigned char vtkType
)
{
    dgEquiSpacedTessellation::VtkSubCell cell;
    cell.nPoints = localVerts.size();
    cell.vtkType = vtkType;

    forAll(cell.pointLabels, pointI)
    {
        cell.pointLabels[pointI] = -1;
    }

    label vertI = 0;

    for (const label localVert : localVerts)
    {
        cell.pointLabels[vertI] = localVert;
        ++vertI;
    }

    cells.append(cell);
}

} // End anonymous namespace


dgEquiSpacedTessellation::dgEquiSpacedTessellation
(
    const label n,
    const label pOrder
)
:
    n_(n),
    pOrder_(pOrder)
{
    if (n_ < 1)
    {
        FatalErrorInFunction
            << "Expected n >= 1, got " << n_
            << exit(FatalError);
    }

    if (pOrder_ < 0)
    {
        FatalErrorInFunction
            << "Expected pOrder >= 0, got " << pOrder_
            << exit(FatalError);
    }
}


void dgEquiSpacedTessellation::computeBasisAt
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
            computeHexBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        case dgCellType::PRISM:
            computePrismBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        case dgCellType::PYRAMID:
            computePyramidBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        case dgCellType::TET:
            computeTetBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        default:
            FatalErrorInFunction
                << "Unsupported DG cell type " << type
                << exit(FatalError);
    }
}


vector dgEquiSpacedTessellation::canonicalEtaFromXi
(
    const vector& xi,
    const dgCellType type
)
{
    const scalar xi2 = xi.y();
    const scalar xi3 = xi.z();

    switch (type)
    {
        case dgCellType::HEX:
            return xiToEtaHex(xi);

        case dgCellType::PRISM:
        {
            const scalar denom = 1.0 - xi3;

            if (mag(denom) <= SMALL)
            {
                return vector(-1.0, xi2, xi3);
            }

            return xiToEtaPrism(xi);
        }

        case dgCellType::PYRAMID:
        {
            const scalar denom = 1.0 - xi3;

            if (mag(denom) <= SMALL)
            {
                return vector(-1.0, -1.0, xi3);
            }

            return xiToEtaPyramid(xi);
        }

        case dgCellType::TET:
        {
            const scalar denom1 = -xi2 - xi3;

            if (mag(denom1) <= SMALL)
            {
                return vector(-1.0, -1.0, xi3);
            }

            return xiToEtaTet(xi);
        }

        default:
            FatalErrorInFunction
                << "Unsupported DG cell type " << type
                << exit(FatalError);
    }

    return vector::zero;
}


dgEquiSpacedTessellation::SamplePoint dgEquiSpacedTessellation::makeSample
(
    const vector& eta,
    const dgCellType type
) const
{
    SamplePoint sample;
    sample.eta = eta;
    computeBasisAt(type, eta, pOrder_, sample.basis);
    return sample;
}


dgCellType dgEquiSpacedTessellation::cellType(const label nPoints)
{
    switch (nPoints)
    {
        case 4:
            return dgCellType::TET;

        case 5:
            return dgCellType::PYRAMID;

        case 6:
            return dgCellType::PRISM;

        case 8:
            return dgCellType::HEX;

        default:
            FatalErrorInFunction
                << "Unsupported DG cell with " << nPoints << " points."
                << exit(FatalError);
    }

    return dgCellType::INVALID;
}


dgEquiSpacedTessellation::CellInfo dgEquiSpacedTessellation::cellInfo
(
    const cellShape& shape,
    const UList<point>& meshPoints
) const
{
    CellInfo info;
    info.vertices = shape.points(meshPoints);
    info.type = cellType(info.vertices.size());

    if (info.type == dgCellType::PRISM)
    {
        info.vertices = orientPrismVertices(info.vertices);
    }

    return info;
}


const List<dgEquiSpacedTessellation::SamplePoint>&
dgEquiSpacedTessellation::samplePoints
(
    const dgCellType type
) const
{
    switch (type)
    {
        case dgCellType::HEX:
            if (!hexSamples_.size())
            {
                buildHexSamples(hexSamples_);
            }
            return hexSamples_;

        case dgCellType::PRISM:
            if (!prismSamples_.size())
            {
                buildPrismSamples(prismSamples_);
            }
            return prismSamples_;

        case dgCellType::PYRAMID:
            if (!pyramidSamples_.size())
            {
                buildPyramidSamples(pyramidSamples_);
            }
            return pyramidSamples_;

        case dgCellType::TET:
            if (!tetSamples_.size())
            {
                buildTetSamples(tetSamples_);
            }
            return tetSamples_;

        default:
            FatalErrorInFunction
                << "Unsupported DG cell type " << type
                << exit(FatalError);
    }

    return hexSamples_;
}


const List<dgEquiSpacedTessellation::SamplePoint>&
dgEquiSpacedTessellation::collapsedSamplePoints
(
    const dgCellType type,
    const label collapseDir
) const
{
    switch (type)
    {
        case dgCellType::HEX:
            if (collapseDir < 0 || collapseDir > 2)
            {
                FatalErrorInFunction
                    << "Expected hexahedron collapseDir in [0,2], got "
                    << collapseDir << exit(FatalError);
            }

            if (!hexCollapsedSamples_[collapseDir].size())
            {
                buildHexCollapsedSamples
                (
                    hexCollapsedSamples_[collapseDir],
                    collapseDir
                );
            }
            return hexCollapsedSamples_[collapseDir];

        case dgCellType::PRISM:
            if (collapseDir != 1)
            {
                FatalErrorInFunction
                    << "Prism collapseDir must be 1 (eta2), got "
                    << collapseDir << exit(FatalError);
            }

            if (!prismCollapsedSamples_.size())
            {
                buildPrismCollapsedSamples(prismCollapsedSamples_);
            }
            return prismCollapsedSamples_;

        default:
            FatalErrorInFunction
                << "Collapsed VTU sampling is only implemented for HEX and "
                << "PRISM cells, "
                << "got cell type " << type << exit(FatalError);
    }

    return hexCollapsedSamples_[0];
}


const List<dgEquiSpacedTessellation::VtkSubCell>&
dgEquiSpacedTessellation::subCells
(
    const dgCellType type
) const
{
    switch (type)
    {
        case dgCellType::HEX:
            if (!hexSubCells_.size())
            {
                buildHexSubCells(hexSubCells_);
            }
            return hexSubCells_;

        case dgCellType::PRISM:
            if (!prismSubCells_.size())
            {
                buildPrismSubCells(prismSubCells_);
            }
            return prismSubCells_;

        case dgCellType::PYRAMID:
            if (!pyramidSubCells_.size())
            {
                buildPyramidSubCells(pyramidSubCells_);
            }
            return pyramidSubCells_;

        case dgCellType::TET:
            if (!tetSubCells_.size())
            {
                buildTetSubCells(tetSubCells_);
            }
            return tetSubCells_;

        default:
            FatalErrorInFunction
                << "Unsupported DG cell type " << type
                << exit(FatalError);
    }

    return hexSubCells_;
}


const List<dgEquiSpacedTessellation::VtkSubCell>&
dgEquiSpacedTessellation::collapsedSubCells
(
    const dgCellType type,
    const label collapseDir
) const
{
    switch (type)
    {
        case dgCellType::HEX:
            if (collapseDir < 0 || collapseDir > 2)
            {
                FatalErrorInFunction
                    << "Expected hexahedron collapseDir in [0,2], got "
                    << collapseDir << exit(FatalError);
            }

            if (!hexCollapsedSubCells_[collapseDir].size())
            {
                buildHexCollapsedSubCells
                (
                    hexCollapsedSubCells_[collapseDir],
                    collapseDir
                );
            }
            return hexCollapsedSubCells_[collapseDir];

        case dgCellType::PRISM:
            if (collapseDir != 1)
            {
                FatalErrorInFunction
                    << "Prism collapseDir must be 1 (eta2), got "
                    << collapseDir << exit(FatalError);
            }

            if (!prismCollapsedSubCells_.size())
            {
                buildPrismCollapsedSubCells(prismCollapsedSubCells_);
            }
            return prismCollapsedSubCells_;

        default:
            FatalErrorInFunction
                << "Collapsed VTU sampling is only implemented for HEX and "
                << "PRISM cells, "
                << "got cell type " << type << exit(FatalError);
    }

    return hexCollapsedSubCells_[0];
}


void dgEquiSpacedTessellation::buildHexSamples(List<SamplePoint>& samples) const
{
    DynamicList<SamplePoint> dynSamples((n_ + 1)*(n_ + 1)*(n_ + 1));

    for (label k = 0; k <= n_; ++k)
    {
        const scalar eta3 = -1.0 + 2.0*scalar(k)/scalar(n_);

        for (label j = 0; j <= n_; ++j)
        {
            const scalar eta2 = -1.0 + 2.0*scalar(j)/scalar(n_);

            for (label i = 0; i <= n_; ++i)
            {
                const scalar eta1 = -1.0 + 2.0*scalar(i)/scalar(n_);
                dynSamples.append(makeSample(vector(eta1, eta2, eta3), dgCellType::HEX));
            }
        }
    }

    samples.transfer(dynSamples);
}


void dgEquiSpacedTessellation::buildPrismSamples(List<SamplePoint>& samples) const
{
    const label nSamples = (n_ + 1)*(n_ + 1)*(n_ + 2)/2;
    DynamicList<SamplePoint> dynSamples(nSamples);

    for (label l = 0; l <= n_; ++l)
    {
        const scalar xi2 = -1.0 + 2.0*scalar(l)/scalar(n_);

        for (label i = 0; i <= n_; ++i)
        {
            const scalar xi3 = -1.0 + 2.0*scalar(i)/scalar(n_);

            for (label j = 0; j <= n_ - i; ++j)
            {
                const scalar xi1 = -1.0 + 2.0*scalar(j)/scalar(n_);
                const vector xi(xi1, xi2, xi3);

                dynSamples.append
                (
                    makeSample
                    (
                        canonicalEtaFromXi(xi, dgCellType::PRISM),
                        dgCellType::PRISM
                    )
                );
            }
        }
    }

    samples.transfer(dynSamples);
}


void dgEquiSpacedTessellation::buildPyramidSamples(List<SamplePoint>& samples) const
{
    const label nSamples = (n_ + 1)*(n_ + 2)*(2*n_ + 3)/6;
    DynamicList<SamplePoint> dynSamples(nSamples);

    for (label k = 0; k <= n_; ++k)
    {
        const scalar xi3 = -1.0 + 2.0*scalar(k)/scalar(n_);
        const label maxIndex = n_ - k;

        for (label j = 0; j <= maxIndex; ++j)
        {
            const scalar xi2 = -1.0 + 2.0*scalar(j)/scalar(n_);

            for (label i = 0; i <= maxIndex; ++i)
            {
                const scalar xi1 = -1.0 + 2.0*scalar(i)/scalar(n_);
                const vector xi(xi1, xi2, xi3);

                dynSamples.append
                (
                    makeSample
                    (
                        canonicalEtaFromXi(xi, dgCellType::PYRAMID),
                        dgCellType::PYRAMID
                    )
                );
            }
        }
    }

    samples.transfer(dynSamples);
}


void dgEquiSpacedTessellation::buildTetSamples(List<SamplePoint>& samples) const
{
    const label nSamples = (n_ + 1)*(n_ + 2)*(n_ + 3)/6;
    DynamicList<SamplePoint> dynSamples(nSamples);

    for (label k = 0; k <= n_; ++k)
    {
        const scalar xi3 = -1.0 + 2.0*scalar(k)/scalar(n_);

        for (label i = 0; i <= n_ - k; ++i)
        {
            const scalar xi2 = -1.0 + 2.0*scalar(i)/scalar(n_);

            for (label j = 0; j <= n_ - k - i; ++j)
            {
                const scalar xi1 = -1.0 + 2.0*scalar(j)/scalar(n_);
                const vector xi(xi1, xi2, xi3);

                dynSamples.append
                (
                    makeSample
                    (
                        canonicalEtaFromXi(xi, dgCellType::TET),
                        dgCellType::TET
                    )
                );
            }
        }
    }

    samples.transfer(dynSamples);
}


void dgEquiSpacedTessellation::buildHexCollapsedSamples
(
    List<SamplePoint>& samples,
    const label collapseDir
) const
{
    DynamicList<SamplePoint> dynSamples((n_ + 1)*(n_ + 1));

    for (label j = 0; j <= n_; ++j)
    {
        const scalar coordJ = -1.0 + 2.0*scalar(j)/scalar(n_);

        for (label i = 0; i <= n_; ++i)
        {
            const scalar coordI = -1.0 + 2.0*scalar(i)/scalar(n_);
            vector eta(vector::zero);

            switch (collapseDir)
            {
                case 0:
                    eta = vector(0.0, coordI, coordJ);
                    break;

                case 1:
                    eta = vector(coordI, 0.0, coordJ);
                    break;

                case 2:
                    eta = vector(coordI, coordJ, 0.0);
                    break;

                default:
                    FatalErrorInFunction
                        << "Expected hexahedron collapseDir in [0,2], got "
                        << collapseDir << exit(FatalError);
            }

            dynSamples.append(makeSample(eta, dgCellType::HEX));
        }
    }

    samples.transfer(dynSamples);
}


void dgEquiSpacedTessellation::buildPrismCollapsedSamples
(
    List<SamplePoint>& samples
) const
{
    const label nSamples = trianglePointCount(n_);
    DynamicList<SamplePoint> dynSamples(nSamples);
    const scalar xi2 = 0.0;

    for (label i = 0; i <= n_; ++i)
    {
        const scalar xi3 = -1.0 + 2.0*scalar(i)/scalar(n_);

        for (label j = 0; j <= n_ - i; ++j)
        {
            const scalar xi1 = -1.0 + 2.0*scalar(j)/scalar(n_);
            const vector xi(xi1, xi2, xi3);

            dynSamples.append
            (
                makeSample
                (
                    canonicalEtaFromXi(xi, dgCellType::PRISM),
                    dgCellType::PRISM
                )
            );
        }
    }

    samples.transfer(dynSamples);
}


void dgEquiSpacedTessellation::buildHexSubCells(List<VtkSubCell>& cells) const
{
    DynamicList<VtkSubCell> dynCells(n_*n_*n_);

    auto hexPoint =
    [&](const label i, const label j, const label k)
    {
        return i + (n_ + 1)*(j + (n_ + 1)*k);
    };

    for (label k = 0; k < n_; ++k)
    {
        for (label j = 0; j < n_; ++j)
        {
            for (label i = 0; i < n_; ++i)
            {
                appendSubCell
                (
                    dynCells,
                    {
                        hexPoint(i, j, k),
                        hexPoint(i + 1, j, k),
                        hexPoint(i + 1, j + 1, k),
                        hexPoint(i, j + 1, k),
                        hexPoint(i, j, k + 1),
                        hexPoint(i + 1, j, k + 1),
                        hexPoint(i + 1, j + 1, k + 1),
                        hexPoint(i, j + 1, k + 1)
                    },
                    vtkHexType
                );
            }
        }
    }

    cells.transfer(dynCells);
}


void dgEquiSpacedTessellation::buildPrismSubCells(List<VtkSubCell>& cells) const
{
    DynamicList<VtkSubCell> dynCells(n_*n_*n_);

    for (label layer = 0; layer < n_; ++layer)
    {
        for (label row = 0; row < n_; ++row)
        {
            for (label i = 0; i < n_ - row; ++i)
            {
                const label A = prismPointIndex(n_, layer, i, row);
                const label B = prismPointIndex(n_, layer, i + 1, row);
                const label C = prismPointIndex(n_, layer, i, row + 1);

                const label a = prismPointIndex(n_, layer + 1, i, row);
                const label b = prismPointIndex(n_, layer + 1, i + 1, row);
                const label c = prismPointIndex(n_, layer + 1, i, row + 1);

                appendSubCell(dynCells, {A, B, C, a, b, c}, vtkWedgeType);

                if (i + row < n_ - 1)
                {
                    const label D = prismPointIndex(n_, layer, i + 1, row + 1);
                    const label d =
                        prismPointIndex(n_, layer + 1, i + 1, row + 1);

                    appendSubCell(dynCells, {B, D, C, b, d, c}, vtkWedgeType);
                }
            }
        }
    }

    cells.transfer(dynCells);
}


void dgEquiSpacedTessellation::buildPyramidSubCells(List<VtkSubCell>& cells) const
{
    DynamicList<VtkSubCell> dynCells(n_*n_*n_);

    for (label layer = 0; layer < n_; ++layer)
    {
        const label m = n_ - layer;

        for (label j = 0; j < m; ++j)
        {
            for (label i = 0; i < m; ++i)
            {
                const label L00 = pyramidPointIndex(n_, layer, i, j);
                const label L10 = pyramidPointIndex(n_, layer, i + 1, j);
                const label L11 = pyramidPointIndex(n_, layer, i + 1, j + 1);
                const label L01 = pyramidPointIndex(n_, layer, i, j + 1);

                const bool hasUpperX = (i < m - 1);
                const bool hasUpperY = (j < m - 1);

                if (hasUpperX && hasUpperY)
                {
                    const label U00 = pyramidPointIndex(n_, layer + 1, i, j);
                    const label U10 = pyramidPointIndex(n_, layer + 1, i + 1, j);
                    const label U11 =
                        pyramidPointIndex(n_, layer + 1, i + 1, j + 1);
                    const label U01 = pyramidPointIndex(n_, layer + 1, i, j + 1);

                    appendSubCell
                    (
                        dynCells,
                        {L00, L10, L11, L01, U00, U10, U11, U01},
                        vtkHexType
                    );
                }
                else if (!hasUpperX && !hasUpperY)
                {
                    const label apex = pyramidPointIndex(n_, layer + 1, i, j);
                    appendSubCell
                    (
                        dynCells,
                        {L00, L10, L11, L01, apex},
                        vtkPyramidType
                    );
                }
                else if (!hasUpperX)
                {
                    const label U0 = pyramidPointIndex(n_, layer + 1, i, j);
                    const label U1 = pyramidPointIndex(n_, layer + 1, i, j + 1);
                    appendSubCell(dynCells, {L00, L10, U0, L01, L11, U1}, vtkWedgeType);
                }
                else
                {
                    const label U0 = pyramidPointIndex(n_, layer + 1, i, j);
                    const label U1 = pyramidPointIndex(n_, layer + 1, i + 1, j);
                    appendSubCell(dynCells, {L00, L01, U0, L10, L11, U1}, vtkWedgeType);
                }
            }
        }
    }

    cells.transfer(dynCells);
}


void dgEquiSpacedTessellation::buildTetSubCells(List<VtkSubCell>& cells) const
{
    DynamicList<VtkSubCell> dynCells(n_*n_*n_);

    for (label layer = 0; layer < n_; ++layer)
    {
        const label m = n_ - layer;

        for (label row = 0; row < m; ++row)
        {
            for (label i = 0; i < m - row; ++i)
            {
                const label A = tetPointIndex(n_, layer, i, row);
                const label B = tetPointIndex(n_, layer, i + 1, row);
                const label C = tetPointIndex(n_, layer, i, row + 1);
                const label a = tetPointIndex(n_, layer + 1, i, row);

                appendSubCell(dynCells, {A, B, C, a}, vtkTetType);

                if (i + row < m - 1)
                {
                    const label D = tetPointIndex(n_, layer, i + 1, row + 1);
                    const label b = tetPointIndex(n_, layer + 1, i + 1, row);
                    const label c = tetPointIndex(n_, layer + 1, i, row + 1);

                    appendSubCell(dynCells, {B, D, C, a, b, c}, vtkWedgeType);
                }
            }
        }
    }

    cells.transfer(dynCells);
}


void dgEquiSpacedTessellation::buildHexCollapsedSubCells
(
    List<VtkSubCell>& cells,
    const label collapseDir
) const
{
    DynamicList<VtkSubCell> dynCells(n_*n_);

    auto quadPoint =
    [&](const label i, const label j)
    {
        return i + (n_ + 1)*j;
    };

    for (label j = 0; j < n_; ++j)
    {
        for (label i = 0; i < n_; ++i)
        {
            const label q00 = quadPoint(i, j);
            const label q10 = quadPoint(i + 1, j);
            const label q11 = quadPoint(i + 1, j + 1);
            const label q01 = quadPoint(i, j + 1);

            if (collapseDir < 0 || collapseDir > 2)
            {
                FatalErrorInFunction
                    << "Expected hexahedron collapseDir in [0,2], got "
                    << collapseDir << exit(FatalError);
            }

            appendSubCell(dynCells, {q00, q10, q11, q01}, vtkQuadType);
        }
    }

    cells.transfer(dynCells);
}


void dgEquiSpacedTessellation::buildPrismCollapsedSubCells
(
    List<VtkSubCell>& cells
) const
{
    DynamicList<VtkSubCell> dynCells(n_*n_);

    for (label row = 0; row < n_; ++row)
    {
        for (label i = 0; i < n_ - row; ++i)
        {
            const label A = triangleRowOffset(n_, row) + i;
            const label B = triangleRowOffset(n_, row) + i + 1;
            const label C = triangleRowOffset(n_, row + 1) + i;

            appendSubCell(dynCells, {A, B, C}, vtkTriType);

            if (i + row < n_ - 1)
            {
                const label D = triangleRowOffset(n_, row + 1) + i + 1;
                appendSubCell(dynCells, {B, D, C}, vtkTriType);
            }
        }
    }

    cells.transfer(dynCells);
}


pointField dgEquiSpacedTessellation::orientPrismVertices
(
    const pointField& rawVertices
)
{
    if (rawVertices.size() != 6)
    {
        FatalErrorInFunction
            << "Expected 6 prism vertices, got " << rawVertices.size()
            << exit(FatalError);
    }

    static const label permutations[6][6] =
    {
        {0, 1, 2, 3, 4, 5},
        {1, 2, 0, 4, 5, 3},
        {2, 0, 1, 5, 3, 4},
        {0, 2, 1, 3, 5, 4},
        {2, 1, 0, 5, 4, 3},
        {1, 0, 2, 4, 3, 5}
    };

    pointField bestVertices(rawVertices);
    scalar bestDet = -GREAT;

    for (label permI = 0; permI < 6; ++permI)
    {
        pointField candidate(6);

        for (label vertI = 0; vertI < 6; ++vertI)
        {
            candidate[vertI] = rawVertices[permutations[permI][vertI]];
        }

        const scalar detJ =
            det(geometricJacobian::PrismJacobian(vector::zero, candidate));

        if (detJ > bestDet)
        {
            bestDet = detJ;
            bestVertices = candidate;
        }
    }

    if (mag(bestDet) <= SMALL)
    {
        FatalErrorInFunction
            << "Degenerate prism encountered while orienting exporter vertices."
            << exit(FatalError);
    }

    if (bestDet < 0)
    {
        WarningInFunction
            << "Prism orientation normalization could not recover a positive "
            << "Jacobian. Exporting with the best available ordering."
            << endl;
    }

    return bestVertices;
}


scalar dgEquiSpacedTessellation::prismOrientationMeasure
(
    const UList<point>& localPoints,
    const FixedList<label, 8>& pointLabels
)
{
    const point& p0 = localPoints[pointLabels[0]];
    const point& p1 = localPoints[pointLabels[1]];
    const point& p2 = localPoints[pointLabels[2]];
    const point& p3 = localPoints[pointLabels[3]];
    const point& p4 = localPoints[pointLabels[4]];
    const point& p5 = localPoints[pointLabels[5]];

    const point baseCentre = (p0 + p1 + p2)/3.0;
    const point topCentre = (p3 + p4 + p5)/3.0;
    const vector baseNormal = (p1 - p0)^(p2 - p0);

    return baseNormal & (topCentre - baseCentre);
}


void dgEquiSpacedTessellation::orientPrismSubCell
(
    const UList<point>& localPoints,
    FixedList<label, 8>& pointLabels,
    const label nPoints
)
{
    if (nPoints != 6)
    {
        return;
    }

    if (prismOrientationMeasure(localPoints, pointLabels) < 0)
    {
        const label tmp12 = pointLabels[1];
        pointLabels[1] = pointLabels[2];
        pointLabels[2] = tmp12;

        const label tmp45 = pointLabels[4];
        pointLabels[4] = pointLabels[5];
        pointLabels[5] = tmp45;
    }
}


void dgEquiSpacedTessellation::appendVtkCells
(
    const dgCellType type,
    const UList<point>& localPoints,
    const label pointOffset,
    DynamicList<label>& connectivity,
    DynamicList<label>& offsets,
    DynamicList<unsigned char>& cellTypes
) const
{
    const List<VtkSubCell>& cells = subCells(type);
    label connectivityOffset = offsets.size() ? offsets[offsets.size() - 1] : 0;

    forAll(cells, cellI)
    {
        FixedList<label, 8> pointLabels = cells[cellI].pointLabels;

        if (type == dgCellType::PRISM)
        {
            orientPrismSubCell(localPoints, pointLabels, cells[cellI].nPoints);
        }

        for (label pointI = 0; pointI < cells[cellI].nPoints; ++pointI)
        {
            connectivity.append(pointOffset + pointLabels[pointI]);
        }

        connectivityOffset += cells[cellI].nPoints;
        offsets.append(connectivityOffset);
        cellTypes.append(cells[cellI].vtkType);
    }
}


void dgEquiSpacedTessellation::appendCollapsedVtkCells
(
    const dgCellType type,
    const label collapseDir,
    const UList<point>& localPoints,
    const label pointOffset,
    DynamicList<label>& connectivity,
    DynamicList<label>& offsets,
    DynamicList<unsigned char>& cellTypes
) const
{
    const List<VtkSubCell>& cells = collapsedSubCells(type, collapseDir);
    label connectivityOffset = offsets.size() ? offsets[offsets.size() - 1] : 0;

    forAll(cells, cellI)
    {
        for (label pointI = 0; pointI < cells[cellI].nPoints; ++pointI)
        {
            connectivity.append(pointOffset + cells[cellI].pointLabels[pointI]);
        }

        connectivityOffset += cells[cellI].nPoints;
        offsets.append(connectivityOffset);
        cellTypes.append(cells[cellI].vtkType);
    }
}


void dgEquiSpacedTessellation::appendVtkCells
(
    const dgCellType type,
    const UList<point>& localPoints,
    const UList<label>& pointMap,
    DynamicList<label>& connectivity,
    DynamicList<label>& offsets,
    DynamicList<unsigned char>& cellTypes
) const
{
    const List<VtkSubCell>& cells = subCells(type);
    label connectivityOffset = offsets.size() ? offsets[offsets.size() - 1] : 0;

    forAll(cells, cellI)
    {
        FixedList<label, 8> pointLabels = cells[cellI].pointLabels;

        if (type == dgCellType::PRISM)
        {
            orientPrismSubCell(localPoints, pointLabels, cells[cellI].nPoints);
        }

        for (label pointI = 0; pointI < cells[cellI].nPoints; ++pointI)
        {
            connectivity.append(pointMap[pointLabels[pointI]]);
        }

        connectivityOffset += cells[cellI].nPoints;
        offsets.append(connectivityOffset);
        cellTypes.append(cells[cellI].vtkType);
    }
}


void dgEquiSpacedTessellation::appendCollapsedVtkCells
(
    const dgCellType type,
    const label collapseDir,
    const UList<point>& localPoints,
    const UList<label>& pointMap,
    DynamicList<label>& connectivity,
    DynamicList<label>& offsets,
    DynamicList<unsigned char>& cellTypes
) const
{
    const List<VtkSubCell>& cells = collapsedSubCells(type, collapseDir);
    label connectivityOffset = offsets.size() ? offsets[offsets.size() - 1] : 0;

    forAll(cells, cellI)
    {
        for (label pointI = 0; pointI < cells[cellI].nPoints; ++pointI)
        {
            connectivity.append(pointMap[cells[cellI].pointLabels[pointI]]);
        }

        connectivityOffset += cells[cellI].nPoints;
        offsets.append(connectivityOffset);
        cellTypes.append(cells[cellI].vtkType);
    }
}

} // End namespace Foam

// ************************************************************************* //
