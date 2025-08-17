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
    A collection of utilities for evaluating orthogonal basis functions and
    their derivatives on reference 2D faces (triangle and quadrilateral),
    for use in Discontinuous Galerkin (DG) methods.

    Supports:
    - Dubiner basis for triangular faces (collapsed coordinate system)
    - Tensor-product Legendre basis for quadrilateral faces

See also

SourceFiles
    basisFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef basisFunctions_H
#define basisFunctions_H

#include "dgCellType.H"
#include "scalar.H"
#include "label.H"
#include "List.H"
#include "error.H"  // For FatalErrorInFunction

namespace Foam
{
namespace math
{

// Return number of basis functions for given cell type and polynomial order
label getNumBasis(const label pOrder, const dgCellType type);

// Evaluate basis functions at a 3D point (eta1, eta2, eta3)
// and their derivatives ∂φ/∂eta1, ∂φ/∂eta2, ∂φ/∂eta3
void computeBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    const dgCellType type,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
);

// Evaluate basis functions and derivatives for hexagonal cells
void computeHexBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
);

// Evaluate basis functions and derivatives for prism cells
void computePrismBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
);

// Evaluate basis functions and derivatives for pyramid cells
void computePyramidBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
);

// Evaluate basis functions and derivatives for tetrahedral cells
void computeTetBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
);

} // namespace math
} // namespace Foam

#endif

// ************************************************************************* //

