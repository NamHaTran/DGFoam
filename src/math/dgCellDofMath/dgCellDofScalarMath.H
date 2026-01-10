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
    dgCellDofScalarMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofScalarMath_H
#define Foam_dgCellDofScalarMath_H

#include "cellDof.H"

namespace Foam
{
// * * * * * * * * * * * * * * * Scalar - Scalar * * * * * * * * * * * * * * * //

inline cellDof<scalar> operator+
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    cellDof<scalar> result(A);
    
    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellDof<scalar> operator-
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    cellDof<scalar> result(A);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellDof<scalar> operator*
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    cellDof<scalar> result(A);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellDof<scalar> operator/
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    cellDof<scalar> result(A);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - Vector * * * * * * * * * * * * * * * //

inline cellDof<vector> operator*
(
    const cellDof<scalar>& A,
    const cellDof<vector>& B
)
{
    cellDof<vector> result(B);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellDof<vector> operator*
(
    const cellDof<vector>& A,
    const cellDof<scalar>& B
)
{
    cellDof<vector> result(A);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellDof<vector> operator/
(
    const cellDof<vector>& A,
    const cellDof<scalar>& B
)
{
    cellDof<vector> result(A);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - tensor * * * * * * * * * * * * * * * //

inline cellDof<tensor> operator*
(
    const cellDof<scalar>& A,
    const cellDof<tensor>& B
)
{
    cellDof<tensor> result(B);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }

    for (label i = 0; i < B.nDof(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellDof<tensor> operator*
(
    const cellDof<tensor>& A,
    const cellDof<scalar>& B
)
{
    return B * A;
}

inline cellDof<tensor> operator/
(
    const cellDof<tensor>& A,
    const cellDof<scalar>& B
)
{
    cellDof<tensor> result(A);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - symmTensor * * * * * * * * * * * * //

inline cellDof<symmTensor> operator*
(
    const cellDof<scalar>& A,
    const cellDof<symmTensor>& B
)
{
    cellDof<symmTensor> result(B);
    
     if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellDof<symmTensor> operator*
(
    const cellDof<symmTensor>& A,
    const cellDof<scalar>& B
)
{
    return B * A;
}

inline cellDof<symmTensor> operator/
(
    const cellDof<symmTensor>& A,
    const cellDof<scalar>& B
)
{
    cellDof<symmTensor> result(A);
    
    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// * * * * * * * * * * * * * * * Scalar - sphericalTensor * * * * * * * * * * //

inline cellDof<sphericalTensor> operator*
(
    const cellDof<scalar>& A,
    const cellDof<sphericalTensor>& B
)
{
    cellDof<sphericalTensor> result(B);
    
    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] * B[i];
    }
    return result;
}

inline cellDof<sphericalTensor> operator*
(
    const cellDof<sphericalTensor>& A,
    const cellDof<scalar>& B
)
{
    return B * A;
}

inline cellDof<sphericalTensor> operator/
(
    const cellDof<sphericalTensor>& A,
    const cellDof<scalar>& B
)
{
    cellDof<sphericalTensor> result(A);

    if (A.nDof() != B.nDof())
    {
        FatalErrorInFunction
            << "Incompatible nDof: " << A.nDof() << " and " << B.nDof()
            << exit(FatalError);
    }

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / B[i];
    }
    return result;
}

// Comparison operators between cellDof<scalar>
inline bool operator<
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    if (A.nDof() != B.nDof()) return false;

    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] < B[i])) return false;
    }
    return true;
}

inline bool operator<=
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    if (A.nDof() != B.nDof()) return false;

    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] <= B[i])) return false;
    }
    return true;
}

inline bool operator>
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    if (A.nDof() != B.nDof()) return false;

    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] > B[i])) return false;
    }
    return true;
}

inline bool operator>=
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    if (A.nDof() != B.nDof()) return false;

    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] >= B[i])) return false;
    }
    return true;
}

inline bool operator==
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    if (A.nDof() != B.nDof()) return false;

    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] == B[i])) return false;
    }
    return true;
}

inline bool operator!=
(
    const cellDof<scalar>& A,
    const cellDof<scalar>& B
)
{
    return !(A == B);
}

} // End namespace Foam

#endif

// ************************************************************************* //
