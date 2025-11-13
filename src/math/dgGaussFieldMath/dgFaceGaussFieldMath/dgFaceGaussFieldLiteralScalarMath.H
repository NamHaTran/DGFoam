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
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{
    
// * * * * * * * * * * * * * * scalar * field * * * * * * * * * * * * * * * //

inline faceGaussField<scalar> operator*
(
    const scalar& s,
    const faceGaussField<scalar>& A
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s * A.plusValue(gp);
        result.minusValueAt(gp) = s * A.minusValue(gp);
    }

    return result;
}

inline faceGaussField<scalar> operator*
(
    const faceGaussField<scalar>& A,
    const scalar& s
)
{
    return s * A;
}

inline faceGaussField<vector> operator*
(
    const scalar& s,
    const faceGaussField<vector>& A
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s * A.plusValue(gp);
        result.minusValueAt(gp) = s * A.minusValue(gp);
    }

    return result;
}

inline faceGaussField<vector> operator*
(
    const faceGaussField<vector>& A,
    const scalar& s
)
{
    return s * A;
}

inline faceGaussField<tensor> operator*
(
    const scalar& s,
    const faceGaussField<tensor>& A
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s * A.plusValue(gp);
        result.minusValueAt(gp) = s * A.minusValue(gp);
    }

    return result;
}

inline faceGaussField<tensor> operator*
(
    const faceGaussField<tensor>& A,
    const scalar& s
)
{
    return s * A;
}

inline faceGaussField<symmTensor> operator*
(
    const scalar& s,
    const faceGaussField<symmTensor>& A
)
{
    faceGaussField<symmTensor> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s * A.plusValue(gp);
        result.minusValueAt(gp) = s * A.minusValue(gp);
    }

    return result;
}

inline faceGaussField<symmTensor> operator*
(
    const faceGaussField<symmTensor>& A,
    const scalar& s
)
{
    return s * A;
}

inline faceGaussField<sphericalTensor> operator*
(
    const scalar& s,
    const faceGaussField<sphericalTensor>& A
)
{
    faceGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s * A.plusValue(gp);
        result.minusValueAt(gp) = s * A.minusValue(gp);
    }

    return result;
}

inline faceGaussField<sphericalTensor> operator*
(
    const faceGaussField<sphericalTensor>& A,
    const scalar& s
)
{
    return s * A;
}

// * * * * * * * * * * * * * * scalar / field * * * * * * * * * * * * * * * //

inline faceGaussField<scalar> operator/
(
    const scalar& s,
    const faceGaussField<scalar>& A
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = s / A.plusValue(gp);
        result.minusValueAt(gp) = s / A.minusValue(gp);
    }

    return result;
}

inline faceGaussField<scalar> operator/
(
    const faceGaussField<scalar>& A,
    const scalar& s
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / s;
        result.minusValueAt(gp) = A.minusValue(gp) / s;
    }

    return result;
}

inline faceGaussField<vector> operator/
(
    const faceGaussField<vector>& A,
    const scalar& s
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / s;
        result.minusValueAt(gp) = A.minusValue(gp) / s;
    }

    return result;
}

inline faceGaussField<tensor> operator/
(
    const faceGaussField<tensor>& A,
    const scalar& s
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / s;
        result.minusValueAt(gp) = A.minusValue(gp) / s;
    }

    return result;
}

inline faceGaussField<symmTensor> operator/
(
    const faceGaussField<symmTensor>& A,
    const scalar& s
)
{
    faceGaussField<symmTensor> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / s;
        result.minusValueAt(gp) = A.minusValue(gp) / s;
    }

    return result;
}

inline faceGaussField<sphericalTensor> operator/
(
    const faceGaussField<sphericalTensor>& A,
    const scalar& s
)
{
    faceGaussField<sphericalTensor> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) / s;
        result.minusValueAt(gp) = A.minusValue(gp) / s;
    }

    return result;
}

// * * * * * * * * * * scalar +- field<scalar> only * * * * * * * * * * * * //

inline faceGaussField<scalar> operator+
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

inline faceGaussField<scalar> operator+
(
    const faceGaussField<scalar>& A,
    const scalar& s
)
{
    return s + A;
}

inline faceGaussField<scalar> operator-
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

inline faceGaussField<scalar> operator-
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

inline faceGaussField<scalar> pow
(
    const faceGaussField<scalar>& A,
    const scalar& s
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = std::pow(A.plusValue(gp),s);
        result.minusValueAt(gp) = std::pow(A.minusValue(gp),s);
    }

    return result;
}

// Comparison operators between faceGaussField<scalar> and literal scalar
inline bool operator>
(
    const faceGaussField<scalar>& A,
    const scalar B
)
{
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        if (!((A.plusValue(gp) > B) && (A.minusValue(gp) > B))) return false;
    }
    return true;
}

inline bool operator>=
(
    const faceGaussField<scalar>& A,
    const scalar B
)
{
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        if (!((A.plusValue(gp) >= B) && (A.minusValue(gp) >= B))) return false;
    }
    return true;
}

inline bool operator<
(
    const faceGaussField<scalar>& A,
    const scalar B
)
{
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        if (!((A.plusValue(gp) < B) && (A.minusValue(gp) < B))) return false;
    }
    return true;
}

inline bool operator<=
(
    const faceGaussField<scalar>& A,
    const scalar B
)
{
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        if (!((A.plusValue(gp) <= B) && (A.minusValue(gp) <= B))) return false;
    }
    return true;
}

inline bool operator==
(
    const faceGaussField<scalar>& A,
    const scalar B
)
{
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        if (!((A.plusValue(gp) == B) && (A.minusValue(gp) == B))) return false;
    }
    return true;
}

inline bool operator!=
(
    const faceGaussField<scalar>& A,
    const scalar B
)
{
    return !(A == B);
}

} // End namespace Foam

#endif

// ************************************************************************* //
