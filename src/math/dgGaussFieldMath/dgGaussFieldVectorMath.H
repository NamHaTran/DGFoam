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
    Operator overloads and utility functions for GaussField<vector>:
    - vector ± vector
    - vector ^ vector
    - mag(vector)
    - vector & tensor / symmTensor / sphericalTensor

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldVectorMath_H
#define Foam_dgGaussFieldVectorMath_H

#include "cellGaussField.H"
#include "GaussField.H"
#include "GaussField.H"
#include "dgFaceGaussFieldVectorMath.H"
#include "dgCellGaussFieldVectorMath.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{

// * * * * * * * * * * * * Vector - Vector * * * * * * * * * * * * //

//- Vector addition
// Component by component
inline GaussField<vector> operator+
(
    const GaussField<vector>& A,
    const GaussField<vector>& B
)
{
    GaussField<vector> result(A);
    result.cellField() = A.cellField() + B.cellField();
    result.faceField() = A.faceField() + B.faceField();
    return result;
}


//- Vector subtraction
// Component by component
inline GaussField<vector> operator-
(
    const GaussField<vector>& A,
    const GaussField<vector>& B
)
{
    GaussField<vector> result(A);
    result.cellField() = A.cellField() - B.cellField();
    result.faceField() = A.faceField() - B.faceField();
    return result;
}

//- Vector cross product
inline GaussField<vector> operator^
(
    const GaussField<vector>& A,
    const GaussField<vector>& B
)
{
    GaussField<vector> result(A);
    result.cellField() = A.cellField() ^ B.cellField();
    result.faceField() = A.faceField() ^ B.faceField();
    return result;
}

//- Vector dot product
inline GaussField<scalar> operator&
(
    const GaussField<vector>& A,
    const GaussField<vector>& B
)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

//- Vector magnitude
inline GaussField<scalar> mag
(
    const GaussField<vector>& A
)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = mag(A.cellField());
    result.faceField() = mag(A.faceField());
    return result;
}

//- Vector squared magnitude
inline GaussField<scalar> magSqr
(
    const GaussField<vector>& A
)
{
    GaussField<scalar> result(A.cellID(), A.dgMesh());
    result.cellField() = magSqr(A.cellField());
    result.faceField() = magSqr(A.faceField());
    return result;
}

// * * * * * * * * * * * * Vector - Tensor * * * * * * * * * * * * //

//- Vector & tensor
// Dot product of vector with each column of the tensor
inline GaussField<vector> operator&
(
    const GaussField<vector>& A,
    const GaussField<tensor>& B
)
{
    GaussField<vector> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

//- Tensor & vector
// Dot product of each row of the tensor with the vector
inline GaussField<vector> operator&
(
    const GaussField<tensor>& A,
    const GaussField<vector>& B
)
{
    GaussField<vector> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}


// * * * * * * * * * * * * Vector - symmTensor * * * * * * * * * * //

//- Vector & symmTensor
// Dot product of vector with each column of the symmTensor
inline GaussField<vector> operator&
(
    const GaussField<vector>& A,
    const GaussField<symmTensor>& B
)
{
    GaussField<vector> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

//- symmTensor & vector
// Dot product of each row of the symmTensor with the vector
inline GaussField<vector> operator&
(
    const GaussField<symmTensor>& A,
    const GaussField<vector>& B
)
{
    GaussField<vector> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}


// * * * * * * * * * * * * Vector - sphericalTensor * * * * * * * * * * //

//- Vector & sphericalTensor
// Dot product of vector with each column of the sphericalTensor
inline GaussField<vector> operator&
(
    const GaussField<vector>& A,
    const GaussField<sphericalTensor>& B
)
{
    GaussField<vector> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

//- sphericalTensor & vector
// Dot product of each row of the sphericalTensor with the vector
inline GaussField<vector> operator&
(
    const GaussField<sphericalTensor>& A,
    const GaussField<vector>& B
)
{
    GaussField<vector> result(A.cellID(), A.dgMesh());
    result.cellField() = A.cellField() & B.cellField();
    result.faceField() = A.faceField() & B.faceField();
    return result;
}

} // End namespace Foam

#endif

// ************************************************************************* //
