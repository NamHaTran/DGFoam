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
    Operator overloads and utility functions for faceGaussField<vector>:
    - vector ± vector
    - vector ^ vector
    - mag(vector)
    - vector & tensor / symmTensor / sphericalTensor

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldVectorMath_H
#define Foam_dgFaceGaussFieldVectorMath_H

#include "faceGaussField.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{
// * * * * * * * * * * Vector - Vector * * * * * * * * * * //
// Sum of 2 vectors
faceGaussField<vector> operator+
(
    const faceGaussField<vector>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) + B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) + B.minusValue(gp);
    }
    return result;
}

// Difference of 2 vectors
faceGaussField<vector> operator-
(
    const faceGaussField<vector>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) - B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) - B.minusValue(gp);
    }
    return result;
}

// Cross product
faceGaussField<vector> operator^
(
    const faceGaussField<vector>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) ^ B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) ^ B.minusValue(gp);
    }
    return result;
}

// Dot product
faceGaussField<scalar> operator&
(
    const faceGaussField<vector>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

// Magnitude

faceGaussField<scalar> mag
(
    const faceGaussField<vector>& A
)
{
    faceGaussField<scalar> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = mag(A.plusValue(gp));
        result.minusValueAt(gp) = mag(A.minusValue(gp));
    }
    return result;
}


// * * * * * * * * * * Vector - Tensor * * * * * * * * * * //

inline faceGaussField<vector> operator&
(
    const faceGaussField<vector>& A,
    const faceGaussField<tensor>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<vector> operator&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}


// * * * * * * * * * * Vector - symmTensor * * * * * * * * * * //

inline faceGaussField<vector> operator&
(
    const faceGaussField<vector>& A,
    const faceGaussField<symmTensor>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<vector> operator&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}


// * * * * * * * * * * Vector - sphericalTensor * * * * * * * * * * //

inline faceGaussField<vector> operator&
(
    const faceGaussField<vector>& A,
    const faceGaussField<sphericalTensor>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

inline faceGaussField<vector> operator&
(
    const faceGaussField<sphericalTensor>& A,
    const faceGaussField<vector>& B
)
{
    faceGaussField<vector> result(A.cellID(), A.dgMesh());
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        result.plusValueAt(gp)  = A.plusValue(gp) & B.plusValue(gp);
        result.minusValueAt(gp) = A.minusValue(gp) & B.minusValue(gp);
    }
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
