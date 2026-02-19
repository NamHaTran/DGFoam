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
    Operator overloads between GaussField<T> types:
    - scalar - scalar
    - scalar - vector
    - scalar - tensor
    - scalar - symmTensor
    - scalar - sphericalTensor

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldScalarMath_H
#define Foam_dgGaussFieldScalarMath_H

#include "cellGaussField.H"
#include "faceGaussField.H"
#include "GaussField.H"
#include "dgFaceGaussFieldScalarMath.H"
#include "dgCellGaussFieldScalarMath.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{
// * * * * * * * * * * Scalar - Scalar * * * * * * * * * * //
inline GaussField<scalar> operator+
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    GaussField<scalar> result(A);
    result.cellField() = A.cellField() + B.cellField();
    result.faceField() = A.faceField() + B.faceField();
    return result;
}

inline GaussField<scalar> operator-
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    GaussField<scalar> result(A);
    result.cellField() = A.cellField() - B.cellField();
    result.faceField() = A.faceField() - B.faceField();
    return result;
}


inline GaussField<scalar> operator*
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    GaussField<scalar> result(A);
    result.cellField() = A.cellField() * B.cellField();
    result.faceField() = A.faceField() * B.faceField();
    return result;
}

inline GaussField<scalar> operator/
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    GaussField<scalar> result(A);
    result.cellField() = A.cellField() / B.cellField();
    result.faceField() = A.faceField() / B.faceField();
    return result;
}


// * * * * * * * * * * Scalar - Vector * * * * * * * * * * //

inline GaussField<vector> operator*
(
    const GaussField<scalar>& A,
    const GaussField<vector>& B
)
{
    GaussField<vector> result(B);
    result.cellField() = A.cellField() * B.cellField();
    result.faceField() = A.faceField() * B.faceField();
    return result;
}

inline GaussField<vector> operator*
(
    const GaussField<vector>& A,
    const GaussField<scalar>& B
)
{
    return B * A;
}

inline GaussField<vector> operator/
(
    const GaussField<vector>& A,
    const GaussField<scalar>& B
)
{
    GaussField<vector> result(A);
    result.cellField() = A.cellField() / B.cellField();
    result.faceField() = A.faceField() / B.faceField();
    return result;
}

// * * * * * * * * * * Scalar - Tensor * * * * * * * * * * //

inline GaussField<tensor> operator*
(
    const GaussField<scalar>& A,
    const GaussField<tensor>& B
)
{
    GaussField<tensor> result(B);
    result.cellField() = A.cellField() * B.cellField();
    result.faceField() = A.faceField() * B.faceField();
    return result;
}

inline GaussField<tensor> operator*
(
    const GaussField<tensor>& A,
    const GaussField<scalar>& B
)
{
    return B * A;
}

inline GaussField<tensor> operator/
(
    const GaussField<tensor>& A,
    const GaussField<scalar>& B
)
{
    GaussField<tensor> result(A);
    result.cellField() = A.cellField() / B.cellField();
    result.faceField() = A.faceField() / B.faceField();
    return result;
}

// * * * * * * * * * * Scalar - symmTensor * * * * * * * * //

inline GaussField<symmTensor> operator*
(
    const GaussField<scalar>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<symmTensor> result(B);
    result.cellField() = A.cellField() * B.cellField();
    result.faceField() = A.faceField() * B.faceField();
    return result;
}

inline GaussField<symmTensor> operator*
(
    const GaussField<symmTensor>& A,
    const GaussField<scalar>& B
)
{
    return B * A;
}

inline GaussField<symmTensor> operator/
(
    const GaussField<symmTensor>& A,
    const GaussField<scalar>& B
)
{
    GaussField<symmTensor> result(A);
    result.cellField() = A.cellField() / B.cellField();
    result.faceField() = A.faceField() / B.faceField();
    return result;
}

// * * * * * * * * * * Scalar - sphericalTensor * * * * * * * * //

inline GaussField<sphericalTensor> operator*
(
    const GaussField<scalar>& A,
    const GaussField<sphericalTensor>& B
)
{
    GaussField<sphericalTensor> result(B);
    result.cellField() = A.cellField() * B.cellField();
    result.faceField() = A.faceField() * B.faceField();
    return result;
}

inline GaussField<sphericalTensor> operator*
(
    const GaussField<sphericalTensor>& A,
    const GaussField<scalar>& B
)
{
    return B * A;
}

inline GaussField<sphericalTensor> operator/
(
    const GaussField<sphericalTensor>& A,
    const GaussField<scalar>& B
)
{
    GaussField<sphericalTensor> result(A);
    result.cellField() = A.cellField() / B.cellField();
    result.faceField() = A.faceField() / B.faceField();
    return result;
}

// Comparison operators between GaussField<scalar>
inline bool operator>
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    return (A.cellField() > B.cellField()) && (A.faceField() > B.faceField());
}

inline bool operator>=
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    return (A.cellField() >= B.cellField()) && (A.faceField() >= B.faceField());
}

inline bool operator<
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    return (A.cellField() < B.cellField()) && (A.faceField() < B.faceField());
}

inline bool operator<=
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    return (A.cellField() <= B.cellField()) && (A.faceField() <= B.faceField());
}

inline bool operator==
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    return (A.cellField() == B.cellField()) && (A.faceField() == B.faceField());
}

inline bool operator!=
(
    const GaussField<scalar>& A,
    const GaussField<scalar>& B
)
{
    return !(A == B);
}

} // End namespace Foam

#endif

// ************************************************************************* //
