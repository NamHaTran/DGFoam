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
    Tensor math operations for boundaryGaussField.

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundaryGaussFieldTensorMath_H
#define Foam_boundaryGaussFieldTensorMath_H

#include "boundaryGaussField.H"
#include "Tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "scalar.H"

namespace Foam
{

inline boundaryGaussField<tensor> dev(const boundaryGaussField<tensor>& A)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = dev(A[i]);
    }
    return result;
}

inline boundaryGaussField<symmTensor> symm(const boundaryGaussField<tensor>& A)
{
    boundaryGaussField<symmTensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = symm(A[i]);
    }
    return result;
}

inline boundaryGaussField<tensor> skew(const boundaryGaussField<tensor>& A)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = skew(A[i]);
    }
    return result;
}

inline boundaryGaussField<scalar> tr(const boundaryGaussField<tensor>& A)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = tr(A[i]);
    }
    return result;
}

inline boundaryGaussField<scalar> det(const boundaryGaussField<tensor>& A)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = det(A[i]);
    }
    return result;
}

inline boundaryGaussField<tensor> T(const boundaryGaussField<tensor>& A)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = T(A[i]);
    }
    return result;
}

inline boundaryGaussField<tensor> inv(const boundaryGaussField<tensor>& A)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = inv(A[i]);
    }
    return result;
}

inline boundaryGaussField<scalar> tr(const boundaryGaussField<symmTensor>& A)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = tr(A[i]);
    }
    return result;
}

inline boundaryGaussField<symmTensor> dev(const boundaryGaussField<symmTensor>& A)
{
    boundaryGaussField<symmTensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = dev(A[i]);
    }
    return result;
}

inline boundaryGaussField<tensor> inv(const boundaryGaussField<symmTensor>& A)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = inv(A[i]);
    }
    return result;
}

inline boundaryGaussField<tensor> operator+
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator-
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator&
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<scalar> operator&&
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] && B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator+
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator-
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator&
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<scalar> operator&&
(
    const boundaryGaussField<tensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] && B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator+
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    return B + A;
}

inline boundaryGaussField<tensor> operator-
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator&
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    return B & A;
}

inline boundaryGaussField<scalar> operator&&
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<tensor>& B
)
{
    return B && A;
}

inline boundaryGaussField<symmTensor> operator+
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<symmTensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline boundaryGaussField<symmTensor> operator-
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<symmTensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline boundaryGaussField<tensor> operator&
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<tensor> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline boundaryGaussField<scalar> operator&&
(
    const boundaryGaussField<symmTensor>& A,
    const boundaryGaussField<symmTensor>& B
)
{
    boundaryGaussField<scalar> result(A.size());
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] && B[i];
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
