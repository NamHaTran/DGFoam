/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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
    Operators between literal scalar and cellDof

SourceFiles
    dgCellDofLiteralScalarMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofLiteralScalarMath_H
#define Foam_dgCellDofLiteralScalarMath_H

#include "cellDof.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include <cmath>

namespace Foam
{
// * * * * * * * * * * * Scalar * Field * * * * * * * * * * * * * * * //

inline cellDof<scalar> operator*
(
    const scalar& a,
    const cellDof<scalar>& A
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<scalar> operator*
(
    const cellDof<scalar>& A,
    const scalar& a
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<vector> operator*
(
    const scalar& a,
    const cellDof<vector>& A
)
{
    cellDof<vector> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<vector> operator*
(
    const cellDof<vector>& A,
    const scalar& a
)
{
    cellDof<vector> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<tensor> operator*
(
    const scalar& a,
    const cellDof<tensor>& A
)
{
    cellDof<tensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<tensor> operator*
(
    const cellDof<tensor>& A,
    const scalar& a
)
{
    cellDof<tensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<symmTensor> operator*
(
    const scalar& a,
    const cellDof<symmTensor>& A
)
{
    cellDof<symmTensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<symmTensor> operator*
(
    const cellDof<symmTensor>& A,
    const scalar& a
)
{
    cellDof<symmTensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<sphericalTensor> operator*
(
    const scalar& a,
    const cellDof<sphericalTensor>& A
)
{
    cellDof<sphericalTensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellDof<sphericalTensor> operator*
(
    const cellDof<sphericalTensor>& A,
    const scalar& a
)
{
    cellDof<sphericalTensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

// * * * * * * * * * * * Scalar / Field * * * * * * * * * * * * * * * //

inline cellDof<scalar> operator/
(
    const scalar& a,
    const cellDof<scalar>& A
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a / A[i];
    }

    return result;
}

inline cellDof<scalar> operator/
(
    const cellDof<scalar>& A,
    const scalar& a
)
{
    cellDof<scalar> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellDof<vector> operator/
(
    const cellDof<vector>& A,
    const scalar& a
)
{
    cellDof<vector> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellDof<tensor> operator/
(
    const cellDof<tensor>& A,
    const scalar& a
)
{
    cellDof<tensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellDof<symmTensor> operator/
(
    const cellDof<symmTensor>& A,
    const scalar& a
)
{
    cellDof<symmTensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellDof<sphericalTensor> operator/
(
    const cellDof<sphericalTensor>& A,
    const scalar& a
)
{
    cellDof<sphericalTensor> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

// * * * * * * * * * Scalar + Field<scalar> * * * * * * * * * * * * * //

inline cellDof<scalar> operator+
(
    const scalar& a,
    const cellDof<scalar>& A
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a + A[i];
    }

    return result;
}

inline cellDof<scalar> operator+
(
    const cellDof<scalar>& A,
    const scalar& a
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a + A[i];
    }

    return result;
}


// * * * * * * * * * Scalar - Field<scalar> * * * * * * * * * * * * * //

inline cellDof<scalar> operator-
(
    const scalar& a,
    const cellDof<scalar>& A
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = a - A[i];
    }

    return result;
}

inline cellDof<scalar> operator-
(
    const cellDof<scalar>& A,
    const scalar& a
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i] - a;
    }

    return result;
}

// pow(field<scalar>, scalar)
inline cellDof<scalar> pow
(
    const cellDof<scalar>& A,
    const scalar& a
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = std::pow(A[i], a);
    }
    return result;
}

// sqrt(field<scalar>)
inline cellDof<scalar> sqrt
(
    const cellDof<scalar>& A
)
{
    cellDof<scalar> result(A);

    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = pow(A[i], 0.5);
    }
    return result;
}

// Comparison operators between literal scalar and cellDof<scalar>
inline bool operator<
(
    const cellDof<scalar>& A,
    const scalar B
)
{
    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] < B)) return false;
    }
    return true;
}

inline bool operator<=
(
    const cellDof<scalar>& A,
    const scalar B
)
{
    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] <= B)) return false;
    }
    return true;
}

inline bool operator>
(
    const cellDof<scalar>& A,
    const scalar B
)
{
    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] > B)) return false;
    }
    return true;
}

inline bool operator>=
(
    const cellDof<scalar>& A,
    const scalar B
)
{
    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] >= B)) return false;
    }
    return true;
}

inline bool operator==
(
    const cellDof<scalar>& A,
    const scalar B
)
{
    for (label i = 0; i < A.nDof(); ++i)
    {
        if (!(A[i] == B)) return false;
    }
    return true;
}

inline bool operator!=
(
    const cellDof<scalar>& A,
    const scalar B
)
{
    return !(A == B);
}

} // End namespace Foam

#endif

// ************************************************************************* //
