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
    dgCellGaussFieldScalarMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldScalarMath_H
#define Foam_dgCellGaussFieldScalarMath_H

#include "cellGaussField.H"

namespace Foam
{
// * * * * * * * * * * * * * * * Scalar - Scalar * * * * * * * * * * * * * * * //

inline cellGaussField<scalar> operator+
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellGaussField<scalar> operator-
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellGaussField<scalar> operator*
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellGaussField<scalar> operator/
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - Vector * * * * * * * * * * * * * * * //

inline cellGaussField<vector> operator*
(
    const cellGaussField<scalar>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellGaussField<vector> operator*
(
    const cellGaussField<vector>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellGaussField<vector> operator/
(
    const cellGaussField<vector>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - tensor * * * * * * * * * * * * * * * //

inline cellGaussField<tensor> operator*
(
    const cellGaussField<scalar>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellGaussField<tensor> operator*
(
    const cellGaussField<tensor>& A,
    const cellGaussField<scalar>& B
)
{
    return B * A;
}

inline cellGaussField<tensor> operator/
(
    const cellGaussField<tensor>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - symmTensor * * * * * * * * * * * * //

inline cellGaussField<symmTensor> operator*
(
    const cellGaussField<scalar>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellGaussField<symmTensor> operator*
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<scalar>& B
)
{
    return B * A;
}

inline cellGaussField<symmTensor> operator/
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - sphericalTensor * * * * * * * * * * //

inline cellGaussField<sphericalTensor> operator*
(
    const cellGaussField<scalar>& A,
    const cellGaussField<sphericalTensor>& B
)
{
    cellGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellGaussField<sphericalTensor> operator*
(
    const cellGaussField<sphericalTensor>& A,
    const cellGaussField<scalar>& B
)
{
    return B * A;
}

inline cellGaussField<sphericalTensor> operator/
(
    const cellGaussField<sphericalTensor>& A,
    const cellGaussField<scalar>& B
)
{
    cellGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// Comparison operators between cellGaussField<scalar>
inline bool operator<
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    if (A.size() != B.size()) return false;

    for (label i = 0; i < A.size(); ++i)
    {
        if (!(A[i] < B[i])) return false;
    }
    return true;
}

inline bool operator<=
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    if (A.size() != B.size()) return false;

    for (label i = 0; i < A.size(); ++i)
    {
        if (!(A[i] <= B[i])) return false;
    }
    return true;
}

inline bool operator>
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    if (A.size() != B.size()) return false;

    for (label i = 0; i < A.size(); ++i)
    {
        if (!(A[i] > B[i])) return false;
    }
    return true;
}

inline bool operator>=
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    if (A.size() != B.size()) return false;

    for (label i = 0; i < A.size(); ++i)
    {
        if (!(A[i] >= B[i])) return false;
    }
    return true;
}

inline bool operator==
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    if (A.size() != B.size()) return false;

    for (label i = 0; i < A.size(); ++i)
    {
        if (!(A[i] == B[i])) return false;
    }
    return true;
}

inline bool operator!=
(
    const cellGaussField<scalar>& A,
    const cellGaussField<scalar>& B
)
{
    return !(A == B);
}

} // End namespace Foam

#endif

// ************************************************************************* //
