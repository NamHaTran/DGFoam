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
