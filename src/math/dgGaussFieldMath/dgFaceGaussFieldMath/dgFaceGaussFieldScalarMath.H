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
// * * * * * * * * * * Scalar - Scalar * * * * * * * * * * //
faceGaussField<scalar> operator+
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) + B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) + B.minusValue(gp);
    }
    return result;
}

faceGaussField<scalar> operator-
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) - B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) - B.minusValue(gp);
    }
    return result;
}


faceGaussField<scalar> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) * B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) * B.minusValue(gp);
    }
    return result;
}

faceGaussField<scalar> operator/
(
    const faceGaussField<scalar>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) / B.minusValue(gp);
    }
    return result;
}


// * * * * * * * * * * Scalar - Vector * * * * * * * * * * //

faceGaussField<vector> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<vector> result(B.cellID(), B.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) * B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) * B.minusValue(gp);
    }
    return result;
}

faceGaussField<vector> operator*
(
    const faceGaussField<vector>& A,
    const faceGaussField<scalar>& B
)
{
    return B * A;
}

faceGaussField<vector> operator/
(
    const faceGaussField<vector>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) / B.minusValue(gp);
    }
    return result;
}


// * * * * * * * * * * Scalar - Tensor * * * * * * * * * * //

faceGaussField<tensor> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<tensor>& B
)
{
    faceGaussField<tensor> result(B.cellID(), B.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) * B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) * B.minusValue(gp);
    }
    return result;
}

faceGaussField<tensor> operator*
(
    const faceGaussField<tensor>& A,
    const faceGaussField<scalar>& B
)
{
    return B * A;
}

faceGaussField<tensor> operator/
(
    const faceGaussField<tensor>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) / B.minusValue(gp);
    }
    return result;
}


// * * * * * * * * * * Scalar - symmTensor * * * * * * * * //

faceGaussField<symmTensor> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<symmTensor> result(B.cellID(), B.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) * B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) * B.minusValue(gp);
    }
    return result;
}

faceGaussField<symmTensor> operator*
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<scalar>& B
)
{
    return B * A;
}

faceGaussField<symmTensor> operator/
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) / B.minusValue(gp);
    }
    return result;
}


// * * * * * * * * * * Scalar - sphericalTensor * * * * * * * * //

faceGaussField<sphericalTensor> operator*
(
    const faceGaussField<scalar>& A,
    const faceGaussField<sphericalTensor>& B
)
{
    faceGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) * B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) * B.minusValue(gp);
    }
    return result;
}

faceGaussField<sphericalTensor> operator*
(
    const faceGaussField<sphericalTensor>& A,
    const faceGaussField<scalar>& B
)
{
    return B * A;
}

faceGaussField<sphericalTensor> operator/
(
    const faceGaussField<sphericalTensor>& A,
    const faceGaussField<scalar>& B
)
{
    faceGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) / B.minusValue(gp);
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
