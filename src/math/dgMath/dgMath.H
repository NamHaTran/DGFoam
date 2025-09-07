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

Description
    Mathematical utilities for DG, including:
        - Gauss integration

SourceFiles
    dgMath.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgMath_H
#define Foam_dgMath_H

#include "List.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

namespace math
{
//- Integrate scalar values with Gauss weights -> scalar
//  values, weights must have the same size.
scalar gaussIntegral
(
    const List<scalar>& values,
    const List<scalar>& weights
);

//- Integrate vector values with Gauss weights -> vector
vector gaussIntegral
(
    const List<vector>& values,
    const List<scalar>& weights
);

//- Integrate tensor values with Gauss weights -> tensor
tensor gaussIntegral
(
    const List<tensor>& values,
    const List<scalar>& weights
);

//- Integrate symmTensor values with Gauss weights -> symmTensor
symmTensor gaussIntegral
(
    const List<symmTensor>& values,
    const List<scalar>& weights
);
};

} // End namespace Foam

#endif

// ************************************************************************* //

