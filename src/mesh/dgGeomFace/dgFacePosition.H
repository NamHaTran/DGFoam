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

    All cell types (tetrahedron, prism, pyramid, hexahedron) are collapsed from
    a reference hexahedron. Each face on a physical cell is mapped to a face on
    the reference cube so that Gauss points and basis functions can be reused.

    Collapse transformation:
      - Prism   : collapse η₁         → E ≈ F, G ≈ H
      - Pyramid : collapse η₁, η₂     → E ≈ F ≈ G ≈ H
      - Tet     : collapse η₁, η₂, η₃ → E ≈ F ≈ G ≈ H, and D ≈ C

    OpenFOAM node indexing for hexahedron:

        H(7) ----------- G(6)
        /|              /|
       / |             / |
     E(4) ----------- F(5)
      |  D(3) -------|-- C(2)     <-- collapse η₃ downward
      | /            | /
      |/             |/
     A(0) -----------B(1)

    Face ordering in OpenFOAM (point labels CCW from outside view):
      Face 0 (left)  : 3 0 4 7 → D A E H → DAEH
      Face 1 (right) : 1 2 6 5 → B C G F → BCGF
      Face 2 (front) : 0 1 5 4 → A B F E → ABFE
      Face 3 (back)  : 2 3 7 6 → C D H G → CDHG
      Face 4 (top)   : 4 5 6 7 → E F G H → EFGH
      Face 5 (bottom): 1 0 3 2 → D C B A → BADC

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
    DAEH = 0,  // (face 0)
    BCGF,      // (face 1)
    ABFE,      // (face 2)
    CDHG,      // (face 3)
    EFGH,      // (face 4)
    BADC,      // (face 5)
    NUM_POSITIONS,
    NONE       // No position (used for boundary faces)
};


// Face enum for tetrahedron (collapsed from cube)
enum class dgFacePositionOnTet
{
    DAE = 0,
    BDE,
    ABE,
    BAD,
    NUM_FACES,
    NONE       // No position (used for boundary faces)
};

// Face enum for pyramid
enum class dgFacePositionOnPyramid
{
    BADC = 0,
    DAE,
    CDE,
    BCE,
    ABE,
    NUM_FACES,
    NONE       // No position (used for boundary faces)
};

// Face enum for prism
enum class dgFacePositionOnPrism
{
    BAC = 0,
    EFG,
    CAEG,
    BCGF,
    ABFE,
    NUM_FACES,
    NONE       // No position (used for boundary faces)
};


// Mapping from tetrahedron face to reference hexahedron face
inline dgFacePosition mapFacePositionFromTet(const dgFacePositionOnTet f)
{
    switch (f)
    {
        case dgFacePositionOnTet::DAE: return dgFacePosition::DAEH;
        case dgFacePositionOnTet::BDE: return dgFacePosition::BCGF;
        case dgFacePositionOnTet::ABE: return dgFacePosition::ABFE;
        case dgFacePositionOnTet::BAD: return dgFacePosition::BADC;
        case dgFacePositionOnTet::NONE: return dgFacePosition::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePositionOnTet enum" << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::DAEH; // Should never reach here
}


// Mapping from pyramid face to reference hexahedron face
inline dgFacePosition mapFacePositionFromPyramid(const dgFacePositionOnPyramid f)
{
    switch (f)
    {
        case dgFacePositionOnPyramid::DAE:  return dgFacePosition::DAEH;
        case dgFacePositionOnPyramid::CDE:  return dgFacePosition::CDHG;
        case dgFacePositionOnPyramid::BCE:  return dgFacePosition::BCGF;
        case dgFacePositionOnPyramid::ABE:  return dgFacePosition::ABFE;
        case dgFacePositionOnPyramid::BADC: return dgFacePosition::BADC;
        case dgFacePositionOnPyramid::NONE: return dgFacePosition::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePositionOnPyramid enum" << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::DAEH; // Should never reach here
}


// Mapping from prism face to reference hexahedron face
inline dgFacePosition mapFacePositionFromPrism(const dgFacePositionOnPrism f)
{
    switch (f)
    {
        case dgFacePositionOnPrism::BAC:   return dgFacePosition::BADC;
        case dgFacePositionOnPrism::EFG:   return dgFacePosition::EFGH;
        case dgFacePositionOnPrism::CAEG:  return dgFacePosition::DAEH;
        case dgFacePositionOnPrism::BCGF:  return dgFacePosition::BCGF;
        case dgFacePositionOnPrism::ABFE:  return dgFacePosition::ABFE;
        case dgFacePositionOnPrism::NONE:  return dgFacePosition::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePositionOnPrism enum" << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::BADC; // Should never reach here
}


// Overload << operator for dgFacePosition
inline Ostream& operator<<(Ostream& os, const dgFacePosition type)
{
    switch (type)
    {
        case dgFacePosition::BADC: return os << "BADC";
        case dgFacePosition::EFGH: return os << "EFGH";
        case dgFacePosition::ABFE: return os << "ABFE";
        case dgFacePosition::BCGF: return os << "BCGF";
        case dgFacePosition::CDHG: return os << "CDHG";
        case dgFacePosition::DAEH: return os << "DAEH";
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
        case 0: return dgFacePositionOnTet::DAE;
        case 1: return dgFacePositionOnTet::BDE;
        case 2: return dgFacePositionOnTet::ABE;
        case 3: return dgFacePositionOnTet::BAD;
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
        case 0: return dgFacePositionOnPrism::BAC;
        case 1: return dgFacePositionOnPrism::EFG;
        case 2: return dgFacePositionOnPrism::CAEG;
        case 3: return dgFacePositionOnPrism::BCGF;
        case 4: return dgFacePositionOnPrism::ABFE;
        default:
            FatalErrorInFunction << "Invalid face index for prism: " << i << abort(FatalError);
    }

    // Dummy
    return dgFacePositionOnPrism::BAC; // Should never reach here
}

// Convert face ID to position on pyramid
inline dgFacePositionOnPyramid convertIdToPositionOnPyramid(label i)
{
    switch (i)
    {
        case 0: return dgFacePositionOnPyramid::BADC;
        case 1: return dgFacePositionOnPyramid::DAE;
        case 2: return dgFacePositionOnPyramid::CDE;
        case 3: return dgFacePositionOnPyramid::BCE;
        case 4: return dgFacePositionOnPyramid::ABE;
        default:
            FatalErrorInFunction << "Invalid face index for pyramid: " << i << abort(FatalError);
    }

    // Dummy
    return dgFacePositionOnPyramid::BADC; // Should never reach here
}

// Convert face ID to position on hexahedron
inline dgFacePosition convertIdToPositionOnHex(label i)
{
    switch (i)
    {
        case 0: return dgFacePosition::DAEH;
        case 1: return dgFacePosition::BCGF;
        case 2: return dgFacePosition::ABFE;
        case 3: return dgFacePosition::CDHG;
        case 4: return dgFacePosition::EFGH;
        case 5: return dgFacePosition::BADC;
        default:
            FatalErrorInFunction << "Invalid face index for hex: " << i << abort(FatalError);
    }

    // Dummy
    return dgFacePosition::DAEH; // Should never reach here
}

} // End namespace Foam

#endif // Foam_dgFacePosition_H

// ************************************************************************* //
