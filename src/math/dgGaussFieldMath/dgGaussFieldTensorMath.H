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
    Tensor math operations for faceGaussField:
    - dev, symm, skew, tr, det, transpose
    - tensor-tensor, tensor-symmTensor, symmTensor-symmTensor operators

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldTensorMath_H
#define Foam_dgFaceGaussFieldTensorMath_H

#include "faceGaussField.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "scalar.H"

namespace Foam
{
namespace dgGaussFieldMath
{

// * * * * * * * * * * * Tensor utility functions * * * * * * * * * * * //

//- Deviatoric part of tensor field
faceGaussField<tensor> dev(const faceGaussField<tensor>& A);

//- Symmetric part of tensor field
faceGaussField<tensor> symm(const faceGaussField<tensor>& A);

//- Skew-symmetric part of tensor field
faceGaussField<tensor> skew(const faceGaussField<tensor>& A);

//- Transpose of tensor field
faceGaussField<tensor> T(const faceGaussField<tensor>& A);

//- Trace of tensor field
faceGaussField<scalar> tr(const faceGaussField<tensor>& A);

//- Determinant of tensor field
faceGaussField<scalar> det(const faceGaussField<tensor>& A);


// * * * * * * * * * * * symmTensor utility functions * * * * * * * * * * * //

//- Trace of symmTensor field
faceGaussField<scalar> tr(const faceGaussField<symmTensor>& A);

//- Deviatoric part of symmTensor field
faceGaussField<symmTensor> dev(const faceGaussField<symmTensor>& A);


// * * * * * * * * * * * tensor-tensor operators * * * * * * * * * * * * * //

//- A + B
faceGaussField<tensor> operator+
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
);

//- A - B
faceGaussField<tensor> operator-
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
);

//- A & B
faceGaussField<tensor> operator&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
);

//- A && B
faceGaussField<scalar> operator&&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<tensor>& B
);

// tensor ± symmTensor
faceGaussField<tensor> operator+
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
);

faceGaussField<tensor> operator-
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
);

faceGaussField<tensor> operator&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
);

faceGaussField<scalar> operator&&
(
    const faceGaussField<tensor>& A,
    const faceGaussField<symmTensor>& B
);

// symmTensor ± tensor
faceGaussField<tensor> operator+
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
);

faceGaussField<tensor> operator-
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
);

faceGaussField<tensor> operator&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
);

faceGaussField<scalar> operator&&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<tensor>& B
);

// symmTensor ± symmTensor
faceGaussField<symmTensor> operator+
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
);

faceGaussField<symmTensor> operator-
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
);

faceGaussField<tensor> operator&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
);

faceGaussField<scalar> operator&&
(
    const faceGaussField<symmTensor>& A,
    const faceGaussField<symmTensor>& B
);

} // End namespace dgGaussFieldMath
} // End namespace Foam

#endif

// ************************************************************************* //
