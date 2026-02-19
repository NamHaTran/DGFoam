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
    dgCellDofTensorMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofTensorMath_H
#define Foam_dgCellDofTensorMath_H

#include "cellDof.H"
#include "Tensor.H"

namespace Foam
{
// * * * * * * * * * * * * * * * Tensor utilities * * * * * * * * * * * * * * //

inline cellDof<tensor> dev(const cellDof<tensor>& A)
{
    cellDof<tensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = dev(A[i]);
    }
    return result;
}

inline cellDof<symmTensor> symm(const cellDof<tensor>& A)
{
    cellDof<symmTensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = symm(A[i]);
    }
    return result;
}

inline cellDof<tensor> skew(const cellDof<tensor>& A)
{
    cellDof<tensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = skew(A[i]);
    }
    return result;
}

inline cellDof<scalar> tr(const cellDof<tensor>& A)
{
    cellDof<scalar> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = tr(A[i]);
    }
    return result;
}

inline cellDof<scalar> det(const cellDof<tensor>& A)
{
    cellDof<scalar> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = det(A[i]);
    }
    return result;
}

inline cellDof<tensor> T(const cellDof<tensor>& A)
{
    cellDof<tensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = T(A[i]);
    }
    return result;
}

inline cellDof<tensor> inv(const cellDof<tensor>& A)
{
    cellDof<tensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = inv(A[i]);
    }
    return result;
}

// symmTensor utilities

inline cellDof<scalar> tr(const cellDof<symmTensor>& A)
{
    cellDof<scalar> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = tr(A[i]);
    }
    return result;
}

inline cellDof<symmTensor> dev(const cellDof<symmTensor>& A)
{
    cellDof<symmTensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = dev(A[i]);
    }
    return result;
}

inline cellDof<symmTensor> T(const cellDof<symmTensor>& A)
{
    cellDof<symmTensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = A[i];
    }
    return result;
}

inline cellDof<tensor> inv(const cellDof<symmTensor>& A)
{
    cellDof<tensor> result(A);
    
    for (label i = 0; i < A.nDof(); ++i)
    {
        result[i] = inv(A[i]);
    }
    return result;
}

// * * * * * * * * * * * * * * Tensor-Tensor Operators * * * * * * * * * * * * //

inline cellDof<tensor> operator+
(
    const cellDof<tensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellDof<tensor> operator-
(
    const cellDof<tensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellDof<tensor> operator&
(
    const cellDof<tensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellDof<scalar> operator&&
(
    const cellDof<tensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] && B[i];
    }
    return result;
}

// * * * * * * * * * * * * Tensor - symmTensor Operators * * * * * * * * * * * * //

inline cellDof<tensor> operator+
(
    const cellDof<tensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellDof<tensor> operator-
(
    const cellDof<tensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellDof<tensor> operator&
(
    const cellDof<tensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellDof<scalar> operator&&
(
    const cellDof<tensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] && B[i];
    }
    return result;
}

// * * * * * * * * * * * * symmTensor - Tensor Operators * * * * * * * * * * * * //

inline cellDof<tensor> operator+
(
    const cellDof<symmTensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellDof<tensor> operator-
(
    const cellDof<symmTensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellDof<tensor> operator&
(
    const cellDof<symmTensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellDof<scalar> operator&&
(
    const cellDof<symmTensor>& A,
    const cellDof<tensor>& B
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
        result[i] = A[i] && B[i];
    }
    return result;
}

// * * * * * * * * * * * * symmTensor - symmTensor Operators * * * * * * * * * * * * //

inline cellDof<symmTensor> operator+
(
    const cellDof<symmTensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] + B[i];
    }
    return result;
}

inline cellDof<symmTensor> operator-
(
    const cellDof<symmTensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] - B[i];
    }
    return result;
}

inline cellDof<tensor> operator&
(
    const cellDof<symmTensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] & B[i];
    }
    return result;
}

inline cellDof<scalar> operator&&
(
    const cellDof<symmTensor>& A,
    const cellDof<symmTensor>& B
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
        result[i] = A[i] && B[i];
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
