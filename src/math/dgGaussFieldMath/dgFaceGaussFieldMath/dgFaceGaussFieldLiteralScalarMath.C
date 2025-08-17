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

\*---------------------------------------------------------------------------*/

#include "dgFaceGaussFieldLiteralScalarMath.H"
#include "faceGaussField.H"

namespace Foam
{
namespace dgGaussFieldMath
{

// * * * * * * * * * * * * * * scalar * field * * * * * * * * * * * * * * * //

template<class Type>
faceGaussField<Type> operator*
(
    const scalar& s,
    const faceGaussField<Type>& A
)
{
    faceGaussField<Type> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s * A.plusValue(gp);
        result.minusValueAt(gp) = s * A.minusValue(gp);
    }

    return result;
}

template<class Type>
faceGaussField<Type> operator*
(
    const faceGaussField<Type>& A,
    const scalar& s
)
{
    return s * A;
}


// * * * * * * * * * * * * * * scalar / field * * * * * * * * * * * * * * * //

template<class Type>
faceGaussField<Type> operator/
(
    const scalar& s,
    const faceGaussField<Type>& A
)
{
    faceGaussField<Type> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s / A.plusValue(gp);
        result.minusValueAt(gp) = s / A.minusValue(gp);
    }

    return result;
}

template<class Type>
faceGaussField<Type> operator/
(
    const faceGaussField<Type>& A,
    const scalar& s
)
{
    faceGaussField<Type> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / s;
        result.minusValueAt(gp) = A.minusValue(gp) / s;
    }

    return result;
}


// * * * * * * * * * * scalar +- field<scalar> only * * * * * * * * * * * * //

faceGaussField<scalar> operator+
(
    const scalar& s,
    const faceGaussField<scalar>& A
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s + A.plusValue(gp);
        result.minusValueAt(gp) = s + A.minusValue(gp);
    }

    return result;
}

faceGaussField<scalar> operator+
(
    const faceGaussField<scalar>& A,
    const scalar& s
)
{
    return s + A;
}

faceGaussField<scalar> operator-
(
    const scalar& s,
    const faceGaussField<scalar>& A
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s - A.plusValue(gp);
        result.minusValueAt(gp) = s - A.minusValue(gp);
    }

    return result;
}

faceGaussField<scalar> operator-
(
    const faceGaussField<scalar>& A,
    const scalar& s
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) - s;
        result.minusValueAt(gp) = A.minusValue(gp) - s;
    }

    return result;
}

} // End namespace dgGaussFieldMath
} // End namespace Foam

// ************************************************************************* //
