/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

Description
    Operator overloads for arithmetic between literal scalar and GaussField

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldLiteralScalarMath_H
#define Foam_dgGaussFieldLiteralScalarMath_H

#include "cellGaussField.H"
#include "faceGaussField.H"
#include "GaussField.H"
#include "dgFaceGaussFieldLiteralScalarMath.H"
#include "dgCellGaussFieldLiteralScalarMath.H"

namespace Foam
{
// * * * * * * * * * * * * * * scalar * field * * * * * * * * * * * * * * * //

inline GaussField<scalar> operator*
(
    const scalar& s,
    const GaussField<scalar>& A
)
{
    GaussField<scalar> result(A);
    result.cellField() = s * A.cellField();
    result.faceField() = s * A.faceField();
    return result;
}

inline GaussField<scalar> operator*
(
    const GaussField<scalar>& A,
    const scalar& s
)
{
    return s * A;
}


// * * * * * * * * * * * * * * scalar / field * * * * * * * * * * * * * * * //

inline GaussField<scalar> operator/
(
    const scalar& s,
    const GaussField<scalar>& A
)
{
    GaussField<scalar> result(A);
    result.cellField() = s / A.cellField();
    result.faceField() = s / A.faceField();
    return result;
}

inline GaussField<scalar> operator/
(
    const GaussField<scalar>& A,
    const scalar& s
)
{
    GaussField<scalar> result(A);
    result.cellField() = A.cellField() / s;
    result.faceField() = A.faceField() / s;
    return result;
}


// * * * * * * * * * * scalar +- field<scalar> only * * * * * * * * * * * * //

inline GaussField<scalar> operator+
(
    const scalar& s,
    const GaussField<scalar>& A
)
{
    GaussField<scalar> result(A);
    result.cellField() = s + A.cellField();
    result.faceField() = s + A.faceField();
    return result;
}

inline GaussField<scalar> operator+
(
    const GaussField<scalar>& A,
    const scalar& s
)
{
    return s + A;
}

inline GaussField<scalar> operator-
(
    const scalar& s,
    const GaussField<scalar>& A
)
{
    GaussField<scalar> result(A);
    result.cellField() = s - A.cellField();
    result.faceField() = s - A.faceField();
    return result;
}

inline GaussField<scalar> operator-
(
    const GaussField<scalar>& A,
    const scalar& s
)
{
    GaussField<scalar> result(A);
    result.cellField() = A.cellField() - s;
    result.faceField() = A.faceField() - s;
    return result;
}

inline GaussField<scalar> pow
(
    const GaussField<scalar>& A,
    const scalar& s
)
{
    GaussField<scalar> result(A);
    result.cellField() = pow(A.cellField(), s);
    result.faceField() = pow(A.faceField(), s);
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
