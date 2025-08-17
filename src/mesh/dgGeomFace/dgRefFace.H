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
    Foam::face

Description
    The dgRefFace class represents a reference face element used in
    the Discontinuous Galerkin (DG) method for 3D finite volume meshes.

    It provides the Gauss integration points (eta1, eta2), weights, and
    polynomial basis functions and their derivatives with respect to eta1
    and eta2 on the reference face.

    This class supports different face types (e.g. triangle, quadrilateral),
    and is initialized based on the polynomial order and face type.

    The computed basis and integration data are later used by the actual
    geometric face classes (dgGeomFace) to evaluate physical fluxes,
    Jacobians, and perform numerical integration over faces.

                            H ----------- G
                           /|          / |
                          / |         /  |
                         E --------- F   |
                         |  |        |   |
                         |  D -------|-- C
                         | /         |  /
                         |/          | /
                         A --------- B

    Standard hex cell node ordering (OpenFOAM cellShape order):
    Vertices:
        0: A    1: B    2: C    3: D
        4: E    5: F    6: G    7: H

    Faces (with standard ordering used in DG and reference coordinates):
        0: A-B-C-D     (bottom)     → eta3 = -1
        1: E-F-G-H     (top)        → eta3 =  1
        2: A-B-F-E     (front)      → eta2 = -1
        3: D-C-G-H     (back)       → eta2 =  1
        4: A-D-H-E     (left)       → eta1 = -1
        5: B-C-G-F     (right)      → eta1 =  1

    Coordinate axis:
         z
         |
         |__ y
        /
       x

    Notes:
    - Each face is written in counter-clockwise order when viewed from outside the cell.
    - This ordering is consistent with standard reference face mappings for tensor-product basis on the reference cube [-1,1]^3.

See also

SourceFiles
    dgRefFace.C

\*---------------------------------------------------------------------------*/

#ifndef dgRefFace_H
#define dgRefFace_H

#include "label.H"
#include "scalar.H"
#include "List.H"
#include "vector.H"
#include "dgFacePosition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgRefFace Declaration
\*---------------------------------------------------------------------------*/

class dgRefFace
{
private:

    // Attributes

        // Polynomial order
        label pOrder_;

        // Number of Gauss points
        label nGauss_;

        // Gauss weights
        List<scalar> wGauss_;

        // Transformed Gauss points for each reference face in reference cube
        List<vector> gaussP_BADC_;
        List<vector> gaussP_EFGH_;
        List<vector> gaussP_ABFE_;
        List<vector> gaussP_CDHG_;
        List<vector> gaussP_BCGF_;
        List<vector> gaussP_DAEH_;

public:

    // Constructors

        //- Construct from polynomial order
        dgRefFace(const label pOrder);

        //- Destructor
        ~dgRefFace() = default;

    // Member Functions

        //- Return polynomial order
        inline label order() const { return pOrder_; }

        //- Return number of Gauss points
        inline label nGauss() const { return nGauss_; }

        //- Return Gauss weights (shared for all faces)
        inline const List<scalar>& weights() const { return wGauss_; }

        //- Return Gauss points for a given face position
        const List<vector>& gaussPoints(const dgFacePosition facePos) const;

        //- Generate Gauss points and weights based on face type and polynomial order
        void generateFaceGaussPointsAndWeights();
};
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
