/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCompressibleZeroGradientBoundaryField_H
#define Foam_dgCompressibleZeroGradientBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class dgCompressibleZeroGradientBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleZeroGradientBoundaryField
 * \brief Zero-gradient compressible DG boundary condition.
 *
 * \details
 * Uses the default base-class ghost-state construction, which simply copies
 * the interior conservative state to the exterior side:
 *
 * \f[
 *     \rho^+ = \rho^-,
 *     \qquad
 *     \rho\mathbf{U}^+ = \rho\mathbf{U}^-,
 *     \qquad
 *     E^+ = E^-.
 * \f]
 *
 * In a finite-volume sense this corresponds to a homogeneous normal-gradient
 * condition on the primitive state.
 */
class dgCompressibleZeroGradientBoundaryField
:
    public dgCompressibleBoundaryField
{
public:

    TypeName("compressible::zeroGradient");

    dgCompressibleZeroGradientBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressibleZeroGradientBoundaryField() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
