/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

Description
    Operator overloads and utility functions for faceGaussField<vector>:
    - vector ± vector
    - vector ^ vector
    - mag(vector)
    - vector & tensor / symmTensor / sphericalTensor

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldVectorMath_H
#define Foam_dgFaceGaussFieldVectorMath_H

#include "faceGaussField.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{
namespace dgGaussFieldMath
{

// * * * * * * * * * * * * Vector - Vector * * * * * * * * * * * * //

//- Vector addition
faceGaussField<vector> operator+
(
    const faceGaussField<vector>& A,
    const faceGaussField<vector>& B
);

//- Vector subtraction
faceGaussField<vector> operator-
(
    const faceGaussField<vector>& A,
    const faceGaussField<vector>& B
);

//- Vector cross product
faceGaussField<vector> operator^
(
    const faceGaussField<vector>& A,
    const faceGaussField<vector>& B
);

//- Vector magnitude
faceGaussField<scalar> mag
(
    const faceGaussField<vector>& A
);


// * * * * * * * * * * * * Vector - Tensor * * * * * * * * * * * * //

//- Vector & tensor
faceGaussField<vector> operator&
(
    const faceGaussField<vector>& A,
    const faceGaussField<tensor>& B
);

//- Tensor & vector
faceGaussField<vector> operator&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<vector>& B
);


// * * * * * * * * * * * * Vector - symmTensor * * * * * * * * * * //

//- Vector & symmTensor
faceGaussField<vector> operator&
(
    const faceGaussField<vector>& A,
    const faceGaussField<symmTensor>& B
);

//- symmTensor & vector
faceGaussField<vector> operator&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<vector>& B
);


// * * * * * * * * * * * * Vector - sphericalTensor * * * * * * * * * * //

//- Vector & sphericalTensor
faceGaussField<vector> operator&
(
    const faceGaussField<vector>& A,
    const faceGaussField<sphericalTensor>& B
);

//- sphericalTensor & vector
faceGaussField<vector> operator&
(
    const faceGaussField<sphericalTensor>& A,
    const faceGaussField<vector>& B
);

} // End namespace dgGaussFieldMath
} // End namespace Foam

#endif

// ************************************************************************* //
