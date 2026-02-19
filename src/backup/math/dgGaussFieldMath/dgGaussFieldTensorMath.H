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
    Tensor math operations for GaussField:
    - dev, symm, skew, tr, det, transpose
    - tensor-tensor, tensor-symmTensor, symmTensor-symmTensor operators

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldTensorMath_H
#define Foam_dgGaussFieldTensorMath_H

#include "GaussField.H"
#include "dgFaceGaussFieldTensorMath.H"
#include "dgCellGaussFieldTensorMath.H"
#include "Tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "scalar.H"

namespace Foam
{
// * * * * * * * * * * * Tensor utility functions * * * * * * * * * * * //

//- Deviatoric part of tensor field
inline GaussField<tensor> dev(const GaussField<tensor>& A)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = dev(A.cellField());
    result.faceField() = dev(A.faceField());
    return result;
}

//- Symmetric part of tensor field
inline GaussField<symmTensor> symm(const GaussField<tensor>& A)
{
    GaussField<symmTensor> result(A.cellID(), A.dgMesh());
    result.cellField() = symm(A.cellField());
    result.faceField() = symm(A.faceField());
    return result;
}

//- Skew-symmetric part of tensor field
inline GaussField<tensor> skew(const GaussField<tensor>& A)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = skew(A.cellField());
    result.faceField() = skew(A.faceField());
    return result;
}

//- Transpose of tensor field
inline GaussField<tensor> T(const GaussField<tensor>& A)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = T(A.cellField());
    result.faceField() = T(A.faceField());
    return result;
}

//- Trace of tensor field
inline GaussField<scalar> tr(const GaussField<tensor>& A)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = tr(A.cellField());
    result.faceField() = tr(A.faceField());
    return result;
}

//- Determinant of tensor field
inline GaussField<scalar> det(const GaussField<tensor>& A)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = det(A.cellField());
    result.faceField() = det(A.faceField());
    return result;
}

//- Inverse of tensor field
inline GaussField<tensor> inv(const GaussField<tensor>& A)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = inv(A.cellField());
    result.faceField() = inv(A.faceField());
    return result;
}

// * * * * * * * * * * * symmTensor utility functions * * * * * * * * * * * //

//- Trace of symmTensor field
inline GaussField<scalar> tr(const GaussField<symmTensor>& A)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = tr(A.cellField());
    result.faceField() = tr(A.faceField());
    return result;
}

//- Deviatoric part of symmTensor field
inline GaussField<symmTensor> dev(const GaussField<symmTensor>& A)
{
    GaussField<symmTensor> result(A.cellID(), A.dgMesh());
    result.cellField() = dev(A.cellField());
    result.faceField() = dev(A.faceField());
    return result;
}

//- Inverse of symmTensor field
inline GaussField<tensor> inv(const GaussField<symmTensor>& A)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = inv(A.cellField());
    result.faceField() = inv(A.faceField());
    return result;
}

// * * * * * * * * * * * tensor-tensor operators * * * * * * * * * * * * * //

//- A + B
inline GaussField<tensor> operator+
(
    const GaussField<tensor>& A,
    const GaussField<tensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() + B.cellField();
    result.faceField() = A.faceField() + B.faceField();
    return result;
}

//- A - B
inline GaussField<tensor> operator-
(
    const GaussField<tensor>& A,
    const GaussField<tensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() - B.cellField();
    result.faceField() = A.faceField() - B.faceField();
    return result;
}

//- A & B
// Normal tensor product: row of A times column of B
inline GaussField<tensor> operator&
(
    const GaussField<tensor>& A,
    const GaussField<tensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

//- A && B
inline GaussField<scalar> operator&&
(
    const GaussField<tensor>& A,
    const GaussField<tensor>& B
)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() && B.cellField();
    result.faceField() = A.faceField() && B.faceField();
    return result;
}

// tensor ± symmTensor
inline GaussField<tensor> operator+
(
    const GaussField<tensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() + B.cellField();
    result.faceField() = A.faceField() + B.faceField();
    return result;
}

inline GaussField<tensor> operator-
(
    const GaussField<tensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() - B.cellField();
    result.faceField() = A.faceField() - B.faceField();
    return result;
}

// Normal tensor product: row of A times column of B
inline GaussField<tensor> operator&
(
    const GaussField<tensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

// Double contraction: element-element product and sum
inline GaussField<scalar> operator&&
(
    const GaussField<tensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() && B.cellField();
    result.faceField() = A.faceField() && B.faceField();
    return result;
}

// symmTensor ± tensor
inline GaussField<tensor> operator+
(
    const GaussField<symmTensor>& A,
    const GaussField<tensor>& B
)
{
    return B + A;
}

inline GaussField<tensor> operator-
(
    const GaussField<symmTensor>& A,
    const GaussField<tensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() - B.cellField();
    result.faceField() = A.faceField() - B.faceField();
    return result;
}

// Normal tensor product: row of A times column of B
inline GaussField<tensor> operator&
(
    const GaussField<symmTensor>& A,
    const GaussField<tensor>& B
)
{
    return B & A;
}

// Double contraction: element-element product and sum
inline GaussField<scalar> operator&&
(
    const GaussField<symmTensor>& A,
    const GaussField<tensor>& B
)
{
    return B && A;
}

// symmTensor ± symmTensor
inline GaussField<symmTensor> operator+
(
    const GaussField<symmTensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<symmTensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() + B.cellField();
    result.faceField() = A.faceField() + B.faceField();
    return result;
}

inline GaussField<symmTensor> operator-
(
    const GaussField<symmTensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<symmTensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() - B.cellField();
    result.faceField() = A.faceField() - B.faceField();
    return result;
}

// Normal tensor product: row of A times column of B
inline GaussField<tensor> operator&
(
    const GaussField<symmTensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<tensor> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

// Double contraction: element-element product and sum
inline GaussField<scalar> operator&&
(
    const GaussField<symmTensor>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() && B.cellField();
    result.faceField() = A.faceField() && B.faceField();
    return result;

}
} // End namespace Foam

#endif

// ************************************************************************* //
