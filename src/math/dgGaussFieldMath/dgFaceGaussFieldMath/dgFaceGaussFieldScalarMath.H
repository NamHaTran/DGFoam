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
    Operator overloads between faceGaussField<T> types:
    - scalar - scalar
    - scalar - vector
    - scalar - tensor
    - scalar - symmTensor
    - scalar - sphericalTensor

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldScalarMath_H
#define Foam_dgFaceGaussFieldScalarMath_H

#include "faceGaussField.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{
namespace dgGaussFieldMath
{

// * * * * * * * * * * Scalar - Scalar * * * * * * * * * * //

//- Add two faceGaussField<scalar>
faceGaussField<scalar> operator+
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
);

//- Subtract two faceGaussField<scalar>
faceGaussField<scalar> operator-
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
);

//- Multiply two faceGaussField<scalar>
faceGaussField<scalar> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
);

//- Divide two faceGaussField<scalar>
faceGaussField<scalar> operator/
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
);


// * * * * * * * * * * Scalar - Vector * * * * * * * * * * //

//- Multiply faceGaussField<scalar> and faceGaussField<vector>
faceGaussField<vector> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<vector>& B
);

//- Multiply faceGaussField<vector> and faceGaussField<scalar>
faceGaussField<vector> operator*
(
    const faceGaussField<vector>& A,
    const faceGaussField<scalar>& B
);


// * * * * * * * * * * Scalar - Tensor * * * * * * * * * * //

//- Multiply faceGaussField<scalar> and faceGaussField<tensor>
faceGaussField<tensor> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<tensor>& B
);

//- Multiply faceGaussField<tensor> and faceGaussField<scalar>
faceGaussField<tensor> operator*
(
    const faceGaussField<tensor>& A,
    const faceGaussField<scalar>& B
);

//- Divide faceGaussField<tensor> by faceGaussField<scalar>
faceGaussField<tensor> operator/
(
    const faceGaussField<tensor>& A,
    const faceGaussField<scalar>& B
);


// * * * * * * * * * * Scalar - symmTensor * * * * * * * * * * //

//- Multiply faceGaussField<scalar> and faceGaussField<symmTensor>
faceGaussField<symmTensor> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<symmTensor>& B
);

//- Multiply faceGaussField<symmTensor> and faceGaussField<scalar>
faceGaussField<symmTensor> operator*
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<scalar>& B
);

//- Divide faceGaussField<symmTensor> by faceGaussField<scalar>
faceGaussField<symmTensor> operator/
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<scalar>& B
);


// * * * * * * * * * * Scalar - sphericalTensor * * * * * * * * * * //

//- Multiply faceGaussField<scalar> and faceGaussField<sphericalTensor>
faceGaussField<sphericalTensor> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<sphericalTensor>& B
);

//- Multiply faceGaussField<sphericalTensor> and faceGaussField<scalar>
faceGaussField<sphericalTensor> operator*
(
    const faceGaussField<sphericalTensor>& A,
    const faceGaussField<scalar>& B
);

//- Divide faceGaussField<sphericalTensor> by faceGaussField<scalar>
faceGaussField<sphericalTensor> operator/
(
    const faceGaussField<sphericalTensor>& A,
    const faceGaussField<scalar>& B
);

} // End namespace dgGaussFieldMath
} // End namespace Foam

#endif

// ************************************************************************* //
