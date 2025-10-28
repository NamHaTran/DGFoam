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
    Operators between literal scalar and cellGaussField

SourceFiles
    dgCellGaussFieldLiteralScalarMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldLiteralScalarMath_H
#define Foam_dgCellGaussFieldLiteralScalarMath_H

#include "cellGaussField.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include <cmath>

namespace Foam
{
// * * * * * * * * * * * Scalar * Field * * * * * * * * * * * * * * * //

inline cellGaussField<scalar> operator*
(
    const scalar& a,
    const cellGaussField<scalar>& A
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<scalar> operator*
(
    const cellGaussField<scalar>& A,
    const scalar& a
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<vector> operator*
(
    const scalar& a,
    const cellGaussField<vector>& A
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<vector> operator*
(
    const cellGaussField<vector>& A,
    const scalar& a
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<tensor> operator*
(
    const scalar& a,
    const cellGaussField<tensor>& A
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<tensor> operator*
(
    const cellGaussField<tensor>& A,
    const scalar& a
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<symmTensor> operator*
(
    const scalar& a,
    const cellGaussField<symmTensor>& A
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<symmTensor> operator*
(
    const cellGaussField<symmTensor>& A,
    const scalar& a
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<sphericalTensor> operator*
(
    const scalar& a,
    const cellGaussField<sphericalTensor>& A
)
{
    cellGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

inline cellGaussField<sphericalTensor> operator*
(
    const cellGaussField<sphericalTensor>& A,
    const scalar& a
)
{
    cellGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a * A[i];
    }

    return result;
}

// * * * * * * * * * * * Scalar / Field * * * * * * * * * * * * * * * //

inline cellGaussField<scalar> operator/
(
    const scalar& a,
    const cellGaussField<scalar>& A
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a / A[i];
    }

    return result;
}

inline cellGaussField<scalar> operator/
(
    const cellGaussField<scalar>& A,
    const scalar& a
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellGaussField<vector> operator/
(
    const scalar& a,
    const cellGaussField<vector>& A
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a / A[i];
    }

    return result;
}

inline cellGaussField<vector> operator/
(
    const cellGaussField<vector>& A,
    const scalar& a
)
{
    cellGaussField<vector> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellGaussField<tensor> operator/
(
    const scalar& a,
    const cellGaussField<tensor>& A
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a / A[i];
    }

    return result;
}

inline cellGaussField<tensor> operator/
(
    const cellGaussField<tensor>& A,
    const scalar& a
)
{
    cellGaussField<tensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellGaussField<symmTensor> operator/
(
    const scalar& a,
    const cellGaussField<symmTensor>& A
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a / A[i];
    }

    return result;
}

inline cellGaussField<symmTensor> operator/
(
    const cellGaussField<symmTensor>& A,
    const scalar& a
)
{
    cellGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

inline cellGaussField<sphericalTensor> operator/
(
    const scalar& a,
    const cellGaussField<sphericalTensor>& A
)
{
    cellGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a / A[i];
    }

    return result;
}

inline cellGaussField<sphericalTensor> operator/
(
    const cellGaussField<sphericalTensor>& A,
    const scalar& a
)
{
    cellGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] / a;
    }

    return result;
}

// * * * * * * * * * Scalar + Field<scalar> * * * * * * * * * * * * * //

inline cellGaussField<scalar> operator+
(
    const scalar& a,
    const cellGaussField<scalar>& A
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a + A[i];
    }

    return result;
}

inline cellGaussField<scalar> operator+
(
    const cellGaussField<scalar>& A,
    const scalar& a
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a + A[i];
    }

    return result;
}


// * * * * * * * * * Scalar - Field<scalar> * * * * * * * * * * * * * //

inline cellGaussField<scalar> operator-
(
    const scalar& a,
    const cellGaussField<scalar>& A
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = a - A[i];
    }

    return result;
}

inline cellGaussField<scalar> operator-
(
    const cellGaussField<scalar>& A,
    const scalar& a
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    

    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = A[i] - a;
    }

    return result;
}

// pow(field<scalar>, scalar)
inline cellGaussField<scalar> pow
(
    const cellGaussField<scalar>& A,
    const scalar& a
)
{
    cellGaussField<scalar> result(A.cellID(), A.dgMesh());
    
    for (label i = 0; i < A.size(); ++i)
    {
        result[i] = std::pow(A[i], a);
    }
    return result;
}
} // End namespace Foam

#endif

// ************************************************************************* //
