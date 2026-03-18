/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2026 Ha Nam Tran
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

Description
    Operator overloads and utility functions for boundaryGaussField<vector>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundaryGaussFieldVectorMath_H
#define Foam_boundaryGaussFieldVectorMath_H

#include "boundaryGaussField.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "VectorSpace.H"

namespace Foam
{

inline boundaryGaussField<vector> operator+
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline boundaryGaussField<vector> operator-
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator*
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline boundaryGaussField<vector> operator^
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] ^ B[i];
    }
    return result;
}

inline boundaryGaussField<scalar> operator&
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<scalar> mag
(
    const boundaryGaussField<vector>& A
)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = mag(A[i]);
    }
    return result;
}

inline boundaryGaussField<scalar> magSqr
(
    const boundaryGaussField<vector>& A
)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = magSqr(A[i]);
    }
    return result;
}

inline boundaryGaussField<vector> operator&
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<tensor>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<vector> operator&
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<vector> operator&
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<vector> operator&
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<vector> operator&
(
    const boundaryGaussField<vector>& A,
    const boundaryGaussField<sphericalTensor>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<vector> operator&
(
    const boundaryGaussField<sphericalTensor>& A,
    const boundaryGaussField<vector>& B
)
{
    boundaryGaussField<vector> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
