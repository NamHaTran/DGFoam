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
    Operator overloads for arithmetic between literal scalar and faceGaussField

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldLiteralScalarMath_H
#define Foam_dgFaceGaussFieldLiteralScalarMath_H

#include "faceGaussField.H"

namespace Foam
{
namespace dgGaussFieldMath
{

// * * * * * * * * * * * * * * scalar * field * * * * * * * * * * * * * * * //

//- Multiply scalar with faceGaussField
template<class Type>
faceGaussField<Type> operator*
(
    const scalar& s,
    const faceGaussField<Type>& A
);

//- Multiply faceGaussField with scalar
template<class Type>
faceGaussField<Type> operator*
(
    const faceGaussField<Type>& A,
    const scalar& s
);


// * * * * * * * * * * * * * * scalar / field * * * * * * * * * * * * * * * //

//- Divide scalar by faceGaussField
template<class Type>
faceGaussField<Type> operator/
(
    const scalar& s,
    const faceGaussField<Type>& A
);

//- Divide faceGaussField by scalar
template<class Type>
faceGaussField<Type> operator/
(
    const faceGaussField<Type>& A,
    const scalar& s
);


// * * * * * * * * * * scalar +- field<scalar> only * * * * * * * * * * * * //

//- Add scalar to faceGaussField<scalar>
faceGaussField<scalar> operator+
(
    const scalar& s,
    const faceGaussField<scalar>& A
);

//- Add faceGaussField<scalar> to scalar
faceGaussField<scalar> operator+
(
    const faceGaussField<scalar>& A,
    const scalar& s
);

//- Subtract faceGaussField<scalar> from scalar
faceGaussField<scalar> operator-
(
    const scalar& s,
    const faceGaussField<scalar>& A
);

//- Subtract scalar from faceGaussField<scalar>
faceGaussField<scalar> operator-
(
    const faceGaussField<scalar>& A,
    const scalar& s
);


} // End namespace dgGaussFieldMath
} // End namespace Foam

#endif

// ************************************************************************* //
