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
    Mathematical operations between cellDof objects

SourceFiles
    dgCellDofVectorMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofVectorMath_H
#define Foam_dgCellDofVectorMath_H

#include "cellDof.H"
#include "VectorSpace.H"

namespace Foam
{
    
// * * * * * * * * * * * * * * * Vector - Vector * * * * * * * * * * * * * * * //

inline cellDof<vector> operator+
(
    const cellDof<vector>& A,
    const cellDof<vector>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellDof<vector> operator-
(
    const cellDof<vector>& A,
    const cellDof<vector>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

// Cross product: A ^ B
inline cellDof<vector> operator^
(
    const cellDof<vector>& A,
    const cellDof<vector>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] ^ B[i];
    }
    return result;
}

// Dot product: A & B
inline cellDof<scalar> operator&
(
    const cellDof<vector>& A,
    const cellDof<vector>& B
)
{
    cellDof<scalar> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

// Magnitude: mag(A)
inline cellDof<scalar> mag
(
    const cellDof<vector>& A
)
{
    cellDof<scalar> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = mag(A[i]);
    }
    return result;
}

// Squared Magnitude: magSqr(A)
inline cellDof<scalar> magSqr
(
    const cellDof<vector>& A
)
{
    cellDof<scalar> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = magSqr(A[i]);
    }
    return result;
}

// Outer product of 2 vectors
inline cellDof<tensor> operator*
(
    const cellDof<vector>& A,
    const cellDof<vector>& B
)
{
    cellDof<tensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Vector - Tensor * * * * * * * * * * * * * * * //

inline cellDof<vector> operator&
(
    const cellDof<vector>& A,
    const cellDof<tensor>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellDof<vector> operator&
(
    const cellDof<tensor>& A,
    const cellDof<vector>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Vector - symmTensor * * * * * * * * * * * * * //

inline cellDof<vector> operator&
(
    const cellDof<vector>& A,
    const cellDof<symmTensor>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellDof<vector> operator&
(
    const cellDof<symmTensor>& A,
    const cellDof<vector>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Vector - sphericalTensor * * * * * * * * * * * //

inline cellDof<vector> operator&
(
    const cellDof<vector>& A,
    const cellDof<sphericalTensor>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellDof<vector> operator&
(
    const cellDof<sphericalTensor>& A,
    const cellDof<vector>& B
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] & B[i];
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
