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

namespace Foam
{
namespace dgGaussFieldMath
{

// Scalar * field
template<class Type>
cellGaussField<Type> operator*
(
    const scalar& a,
    const cellGaussField<Type>& A
);

template<class Type>
cellGaussField<Type> operator*
(
    const cellGaussField<Type>& A,
    const scalar& a
);

// Scalar / field
template<class Type>
cellGaussField<Type> operator/
(
    const scalar& a,
    const cellGaussField<Type>& A
);

template<class Type>
cellGaussField<Type> operator/
(
    const cellGaussField<Type>& A,
    const scalar& a
);

// Scalar + field<scalar>
template<class Type>
cellGaussField<Type> operator+
(
    const scalar& a,
    const cellGaussField<scalar>& A
);

template<class Type>
cellGaussField<Type> operator+
(
    const cellGaussField<scalar>& A,
    const scalar& a
);

// Scalar - field<scalar>
template<class Type>
cellGaussField<Type> operator-
(
    const scalar& a,
    const cellGaussField<scalar>& A
);

template<class Type>
cellGaussField<Type> operator-
(
    const cellGaussField<scalar>& A,
    const scalar& a
);

} // End namespace dgGaussFieldMath
} // End namespace Foam

#endif

// ************************************************************************* //
