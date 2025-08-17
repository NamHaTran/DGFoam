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

#include "dgFaceGaussFieldTensorMath.H"
#include "faceGaussField.H"

namespace Foam
{
namespace dgGaussFieldMath
{

// * * * * * * * * * * * Tensor utilities * * * * * * * * * * //

inline faceGaussField<tensor> dev(const faceGaussField<tensor>& A)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = dev(A.plusValue(gp));
        result.minusValueAt(gp) = dev(A.minusValue(gp));
    }
    return result;
}

inline faceGaussField<tensor> symm(const faceGaussField<tensor>& A)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = symm(A.plusValue(gp));
        result.minusValueAt(gp) = symm(A.minusValue(gp));
    }
    return result;
}

inline faceGaussField<tensor> skew(const faceGaussField<tensor>& A)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = skew(A.plusValue(gp));
        result.minusValueAt(gp) = skew(A.minusValue(gp));
    }
    return result;
}

inline faceGaussField<scalar> tr(const faceGaussField<tensor>& A)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = tr(A.plusValue(gp));
        result.minusValueAt(gp) = tr(A.minusValue(gp));
    }
    return result;
}

inline faceGaussField<scalar> det(const faceGaussField<tensor>& A)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = det(A.plusValue(gp));
        result.minusValueAt(gp) = det(A.minusValue(gp));
    }
    return result;
}

inline faceGaussField<tensor> T(const faceGaussField<tensor>& A)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = T(A.plusValue(gp));
        result.minusValueAt(gp) = T(A.minusValue(gp));
    }
    return result;
}


// * * * * * * * symmTensor utilities * * * * * * * //

inline faceGaussField<scalar> tr(const faceGaussField<symmTensor>& A)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = tr(A.plusValue(gp));
        result.minusValueAt(gp) = tr(A.minusValue(gp));
    }
    return result;
}

inline faceGaussField<symmTensor> dev(const faceGaussField<symmTensor>& A)
{
    faceGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = dev(A.plusValue(gp));
        result.minusValueAt(gp) = dev(A.minusValue(gp));
    }
    return result;
}

// * * * * * * * tensor - tensor * * * * * * * //

inline faceGaussField<tensor> operator+
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) + B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) + B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<tensor> operator-
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) - B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) - B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<tensor> operator&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<scalar> operator&&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) && B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) && B.minusValue(gp);
    }
    return result;
}

// * * * * * * * tensor ± symmTensor * * * * * * * //

inline faceGaussField<tensor> operator+
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) + B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) + B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<tensor> operator-
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) - B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) - B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<tensor> operator&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<scalar> operator&&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) && B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) && B.minusValue(gp);
    }
    return result;
}

// * * * * * * * symmTensor ± tensor * * * * * * * //

inline faceGaussField<tensor> operator+
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
)
{
    return B + A;
}

inline faceGaussField<tensor> operator-
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) - B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) - B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<tensor> operator&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
)
{
    return B & A;
}

inline faceGaussField<scalar> operator&&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
)
{
    return B && A;
}

// * * * * * * * symmTensor ± symmTensor * * * * * * * //

inline faceGaussField<symmTensor> operator+
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) + B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) + B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<symmTensor> operator-
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<symmTensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) - B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) - B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<tensor> operator&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<tensor> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<scalar> operator&&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) && B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) && B.minusValue(gp);
    }
    return result;
}

} // End namespace dgGaussFieldMath
} // End namespace Foam

// ************************************************************************* //
