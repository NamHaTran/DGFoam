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
    Foam::dgGaussFieldMath

Description
    Mathematical operations between cellGaussField objects

SourceFiles
    dgCellGaussFieldTensorMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldTensorMath_H
#define Foam_dgCellGaussFieldTensorMath_H

#include "cellGaussField.H"
#include "Tensor.H"

namespace Foam
{
// * * * * * * * * * * * * * * * Tensor utilities * * * * * * * * * * * * * * //

inline cellGaussField<tensor> dev(const cellGaussField<tensor>& A)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = dev(A[i]);
    }
    return result;
}

inline cellGaussField<symmTensor> symm(const cellGaussField<tensor>& A)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = symm(A[i]);
    }
    return result;
}

inline cellGaussField<tensor> skew(const cellGaussField<tensor>& A)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = skew(A[i]);
    }
    return result;
}

inline cellGaussField<scalar> tr(const cellGaussField<tensor>& A)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = tr(A[i]);
    }
    return result;
}

inline cellGaussField<scalar> det(const cellGaussField<tensor>& A)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = det(A[i]);
    }
    return result;
}

inline cellGaussField<tensor> T(const cellGaussField<tensor>& A)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = T(A[i]);
    }
    return result;
}

inline cellGaussField<tensor> inv(const cellGaussField<tensor>& A)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = inv(A[i]);
    }
    return result;
}

// symmTensor utilities

inline cellGaussField<scalar> tr(const cellGaussField<symmTensor>& A)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = tr(A[i]);
    }
    return result;
}

inline cellGaussField<symmTensor> dev(const cellGaussField<symmTensor>& A)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = dev(A[i]);
    }
    return result;
}

inline cellGaussField<symmTensor> T(const cellGaussField<symmTensor>& A)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i];
    }
    return result;
}

inline cellGaussField<tensor> inv(const cellGaussField<symmTensor>& A)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = inv(A[i]);
    }
    return result;
}

// * * * * * * * * * * * * * * Tensor-Tensor Operators * * * * * * * * * * * * //

inline cellGaussField<tensor> operator+
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator-
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellGaussField<scalar> operator&&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] && B[i];
    }
    return result;
}

// * * * * * * * * * * * * Tensor - symmTensor Operators * * * * * * * * * * * * //

inline cellGaussField<tensor> operator+
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator-
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellGaussField<scalar> operator&&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] && B[i];
    }
    return result;
}

// * * * * * * * * * * * * symmTensor - Tensor Operators * * * * * * * * * * * * //

inline cellGaussField<tensor> operator+
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator-
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellGaussField<scalar> operator&&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] && B[i];
    }
    return result;
}

// * * * * * * * * * * * * symmTensor - symmTensor Operators * * * * * * * * * * * * //

inline cellGaussField<symmTensor> operator+
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellGaussField<symmTensor> operator-
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellGaussField<scalar> operator&&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] && B[i];
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
