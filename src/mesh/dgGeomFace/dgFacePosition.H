/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgFacePosition

Description
    Enumeration for face positions used in Discontinuous Galerkin (DG) methods.

    The face positions are defined based on the vertex ordering convention
    presented in the book:

        Sherwin, S.J., Karniadakis, G.E.
        "Spectral/hp Element Methods for Computational Fluid Dynamics",
        2nd Edition, Oxford University Press, 2005.

    In this reference, each standard element (hexahedron, prism, pyramid,
    tetrahedron) is defined in a reference coordinate system (ξ₁, ξ₂, ξ₃),
    and the faces are identified according to the counter-clockwise (CCW)
    ordering of their vertices when viewed from outside the element.

    The enumeration dgFacePosition encodes the relative position of a face
    within a cell by matching the local vertex ordering of OpenFOAM cells
    to the corresponding face definitions in the reference element used
    in DG formulations.

    The vertex labels (A, B, C, D, E, F, G, H) and the local reference
    coordinate directions (ξ₁, ξ₂, ξ₃) follow exactly the notation used
    in the above reference, as illustrated in the standard-element figures
    for:
        - (a) Hexahedron
        - (b) Prism
        - (c) Pyramid
        - (d) Tetrahedron

    This consistent definition ensures that:
        - Face-based basis functions are evaluated on the correct reference face
        - Face normals and orientation are handled consistently
        - Surface integrals in DG formulations are computed correctly
    

SourceFiles
    dgFacePosition.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFacePosition_H
#define Foam_dgFacePosition_H

#include "error.H"
#include "Ostream.H"

namespace Foam
{

// Reference hexahedron face positions (OpenFOAM CCW ordering)
enum class dgFacePosition
{
    ACEG = 0,  // (face 0)
    BDFH,      // (face 1)
    ABEF,      // (face 2)
    CDGH,      // (face 3)
    ABCD,      // (face 4)
    EFGH,      // (face 5)
    NUM_POSITIONS,
    NONE       // No position (used for boundary faces)
};


// Face enum for tetrahedron (collapsed from cube)
enum class dgFacePositionOnTet
{
    BCE = 0,
    ACE,
    ABE,
    ABC,
    NUM_FACES,
    NONE       // No position (used for boundary faces)
};

// Face enum for pyramid
enum class dgFacePositionOnPyramid
{
    ABCD = 0,
    ACE,
    CDE,
    BDE,
    ABE,
    NUM_FACES,
    NONE       // No position (used for boundary faces)
};

// Face enum for prism
enum class dgFacePositionOnPrism
{
    ABE = 0,
    CDF,
    BCEF,
    ABCD,
    ADEF,
    NUM_FACES,
    NONE       // No position (used for boundary faces)
};


// Mapping from tetrahedron face to reference hexahedron face
inline dgFacePosition mapFacePositionFromTet(const dgFacePositionOnTet f)
{
    switch (f)
    {
        case dgFacePositionOnTet::ACE: return dgFacePosition::ACEG;
        case dgFacePositionOnTet::BCE: return dgFacePosition::BDFH;
        case dgFacePositionOnTet::ABE: return dgFacePosition::ABEF;
        case dgFacePositionOnTet::ABC: return dgFacePosition::ABCD;
        case dgFacePositionOnTet::NONE: return dgFacePosition::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePositionOnTet enum" << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::ABCD; // Should never reach here
}

// Mapping from reference hexahedron face to tetrahedron face
inline dgFacePositionOnTet mapFacePositionToTet(const dgFacePosition f)
{
    switch (f)
    {
        case dgFacePosition::ACEG: return dgFacePositionOnTet::ACE;
        case dgFacePosition::BDFH: return dgFacePositionOnTet::BCE;
        case dgFacePosition::ABEF: return dgFacePositionOnTet::ABE;
        case dgFacePosition::ABCD: return dgFacePositionOnTet::ABC;
        case dgFacePosition::NONE: return dgFacePositionOnTet::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePosition enum" << abort(FatalError);
    } 
    
    // Dummy
    return dgFacePositionOnTet::ABC; // Should never reach here
}


// Mapping from pyramid face to reference hexahedron face
inline dgFacePosition mapFacePositionFromPyramid(const dgFacePositionOnPyramid f)
{
    switch (f)
    {
        case dgFacePositionOnPyramid::ACE:  return dgFacePosition::ACEG;
        case dgFacePositionOnPyramid::CDE:  return dgFacePosition::CDGH;
        case dgFacePositionOnPyramid::BDE:  return dgFacePosition::BDFH;
        case dgFacePositionOnPyramid::ABE:  return dgFacePosition::ABEF;
        case dgFacePositionOnPyramid::ABCD: return dgFacePosition::ABCD;
        case dgFacePositionOnPyramid::NONE: return dgFacePosition::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePositionOnPyramid enum" << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::ABCD; // Should never reach here
}

// Mapping from reference hexahedron face to pyramid face
inline dgFacePositionOnPyramid mapFacePositionToPyramid(const dgFacePosition f)
{
    switch (f)
    {
        case dgFacePosition::ACEG: return dgFacePositionOnPyramid::ACE;
        case dgFacePosition::CDGH: return dgFacePositionOnPyramid::CDE;
        case dgFacePosition::BDFH: return dgFacePositionOnPyramid::BDE;
        case dgFacePosition::ABEF: return dgFacePositionOnPyramid::ABE;
        case dgFacePosition::ABCD: return dgFacePositionOnPyramid::ABCD;
        case dgFacePosition::NONE: return dgFacePositionOnPyramid::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePosition enum" << abort(FatalError);
    }
    // Dummy
    return dgFacePositionOnPyramid::ABCD; // Should never reach here
}


// Mapping from prism face to reference hexahedron face
inline dgFacePosition mapFacePositionFromPrism(const dgFacePositionOnPrism f)
{
    switch (f)
    {
        case dgFacePositionOnPrism::ABE:   return dgFacePosition::ABEF;
        case dgFacePositionOnPrism::CDF:   return dgFacePosition::CDGH;
        case dgFacePositionOnPrism::BCEF:  return dgFacePosition::BDFH;
        case dgFacePositionOnPrism::ABCD:  return dgFacePosition::ABCD;
        case dgFacePositionOnPrism::ADEF:  return dgFacePosition::ADEH;
        case dgFacePositionOnPrism::NONE:  return dgFacePosition::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePositionOnPrism enum" << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::ABCD; // Should never reach here
}

// Mapping from reference hexahedron face to prism face
inline dgFacePositionOnPrism mapFacePositionToPrism(const dgFacePosition f)
{
    switch (f)
    {
        case dgFacePosition::ABEF:  return dgFacePositionOnPrism::ABE;
        case dgFacePosition::CDGH:  return dgFacePositionOnPrism::CDF;
        case dgFacePosition::BDFH:  return dgFacePositionOnPrism::BCEF;
        case dgFacePosition::ABCD:  return dgFacePositionOnPrism::ABCD;
        case dgFacePosition::ADEH:  return dgFacePositionOnPrism::ADEF;
        case dgFacePosition::NONE:  return dgFacePositionOnPrism::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePosition enum" << abort(FatalError);
    }
    // Dummy
    return dgFacePositionOnPrism::ABCD; // Should never reach here
}

// Overload << operator for dgFacePosition
inline Ostream& operator<<(Ostream& os, const dgFacePosition type)
{
    switch (type)
    {
        case dgFacePosition::ACEG: return os << "ACEG";
        case dgFacePosition::BDFH: return os << "BDFH";
        case dgFacePosition::ABEF: return os << "ABEF";
        case dgFacePosition::CDGH: return os << "CDGH";
        case dgFacePosition::ABCD: return os << "ABCD";
        case dgFacePosition::EFGH: return os << "EFGH";
        case dgFacePosition::NONE: return os << "NONE";
        default:                   return os << "UNKNOWN";
    }

    // Dummy
    return os << "UNKNOWN"; // Should never reach here
}

// Convert face ID to position on tetrahedron
inline dgFacePositionOnTet convertIdToPositionOnTet(label i)
{
    switch (i)
    {
        case 0: return dgFacePositionOnTet::BCE;
        case 1: return dgFacePositionOnTet::ACE;
        case 2: return dgFacePositionOnTet::ABE;
        case 3: return dgFacePositionOnTet::ABC;
        default:
            FatalErrorInFunction << "Invalid face index for tet: " << i << abort(FatalError);
    }

    // Dummy
    return dgFacePositionOnTet::DAE; // Should never reach here
}

// Convert face ID to position on pyramid
inline dgFacePositionOnPrism convertIdToPositionOnPrism(label i)
{
    switch (i)
    {
        case 0: return dgFacePositionOnPrism::ABE;
        case 1: return dgFacePositionOnPrism::CDF;
        case 2: return dgFacePositionOnPrism::BCEF;
        case 3: return dgFacePositionOnPrism::ABCD;
        case 4: return dgFacePositionOnPrism::ADEF;
        default:
            FatalErrorInFunction << "Invalid face index for prism: " << i << abort(FatalError);
    }

    // Dummy
    return dgFacePositionOnPrism::ABE; // Should never reach here
}

// Convert face ID to position on pyramid
inline dgFacePositionOnPyramid convertIdToPositionOnPyramid(label i)
{
    switch (i)
    {
        case 0: return dgFacePositionOnPyramid::ABCD;
        case 1: return dgFacePositionOnPyramid::ACE;
        case 2: return dgFacePositionOnPyramid::CDE;
        case 3: return dgFacePositionOnPyramid::BDE;
        case 4: return dgFacePositionOnPyramid::ABE;
        default:
            FatalErrorInFunction << "Invalid face index for pyramid: " << i << abort(FatalError);
    }

    // Dummy
    return dgFacePositionOnPyramid::ABCD; // Should never reach here
}

// Convert face ID to position on hexahedron
inline dgFacePosition convertIdToPositionOnHex(label i)
{
    switch (i)
    {
        case 0: return dgFacePosition::ACEG;
        case 1: return dgFacePosition::BDFH;
        case 2: return dgFacePosition::ABEF;
        case 3: return dgFacePosition::CDGH;
        case 4: return dgFacePosition::ABCD;
        case 5: return dgFacePosition::EFGH;
        default:
            FatalErrorInFunction << "Invalid face index for hex: " << i << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::ABCD; // Should never reach here
}

} // End namespace Foam

#endif // Foam_dgFacePosition_H

// ************************************************************************* //
