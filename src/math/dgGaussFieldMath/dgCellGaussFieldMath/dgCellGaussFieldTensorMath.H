/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Mathematical operations between cellGaussField objects

SourceFiles
    dgCellGaussFieldTensorMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldTensorMath_H
#define Foam_dgCellGaussFieldTensorMath_H

#include "cellGaussField.H"

namespace Foam
{
namespace dgGaussFieldMath
{

// Tensor utilities
cellGaussField<tensor> dev(const cellGaussField<tensor>& A);
cellGaussField<tensor> symm(const cellGaussField<tensor>& A);
cellGaussField<tensor> skew(const cellGaussField<tensor>& A);
cellGaussField<scalar> tr(const cellGaussField<tensor>& A);
cellGaussField<scalar> det(const cellGaussField<tensor>& A);
cellGaussField<tensor> T(const cellGaussField<tensor>& A);

// symmTensor utilities
cellGaussField<scalar> tr(const cellGaussField<symmTensor>& A);
cellGaussField<symmTensor> dev(const cellGaussField<symmTensor>& A);
cellGaussField<symmTensor> T(const cellGaussField<symmTensor>& A);

// Tensor-tensor operators
cellGaussField<tensor> operator+
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
);

cellGaussField<tensor> operator-
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
);

cellGaussField<tensor> operator&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
);

cellGaussField<scalar> operator&&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<tensor>& B
);

// Tensor - symmTensor
cellGaussField<tensor> operator+
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
);

cellGaussField<tensor> operator-
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
);

cellGaussField<tensor> operator&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
);

cellGaussField<scalar> operator&&
(
    const cellGaussField<tensor>& A,
    const cellGaussField<symmTensor>& B
);

// symmTensor - Tensor Operators

inline cellGaussField<tensor> operator+
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
);

inline cellGaussField<tensor> operator-
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
);

inline cellGaussField<tensor> operator&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
);

inline cellGaussField<scalar> operator&&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<tensor>& B
);

// symmTensor - symmTensor
cellGaussField<symmTensor> operator+
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
);

cellGaussField<symmTensor> operator-
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
);

cellGaussField<tensor> operator&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
);

cellGaussField<scalar> operator&&
(
    const cellGaussField<symmTensor>& A,
    const cellGaussField<symmTensor>& B
);

} // End namespace dgGaussFieldMath
} // End namespace Foam

#endif

// ************************************************************************* //
