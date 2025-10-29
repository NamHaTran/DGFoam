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
    dgCellGaussFieldVectorMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldVectorMath_H
#define Foam_dgCellGaussFieldVectorMath_H

#include "cellGaussField.H"

namespace Foam
{
    
// * * * * * * * * * * * * * * * Vector - Vector * * * * * * * * * * * * * * * //

inline cellGaussField<vector> operator+
(
    const cellGaussField<vector>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellGaussField<vector> operator-
(
    const cellGaussField<vector>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

// Cross product: A ^ B
inline cellGaussField<vector> operator^
(
    const cellGaussField<vector>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] ^ B[i];
    }
    return result;
}

// Dot product: A & B
inline cellGaussField<scalar> operator&
(
    const cellGaussField<vector>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

// Magnitude: mag(A)
inline cellGaussField<scalar> mag
(
    const cellGaussField<vector>& A
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = mag(A[i]);
    }
    return result;
}

// * * * * * * * * * * * * * * * Vector - Tensor * * * * * * * * * * * * * * * //

inline cellGaussField<vector> operator&
(
    const cellGaussField<vector>& A,
    const cellGaussField<tensor>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellGaussField<vector> operator&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Vector - symmTensor * * * * * * * * * * * * * //

inline cellGaussField<vector> operator&
(
    const cellGaussField<vector>& A,
    const cellGaussField<symmTensor>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellGaussField<vector> operator&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Vector - sphericalTensor * * * * * * * * * * * //

inline cellGaussField<vector> operator&
(
    const cellGaussField<vector>& A,
    const cellGaussField<sphericalTensor>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellGaussField<vector> operator&
(
    const cellGaussField<sphericalTensor>& A,
    const cellGaussField<vector>& B
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
