/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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

#ifndef Foam_dgCompressibleSlipWallBoundaryField_H
#define Foam_dgCompressibleSlipWallBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dgCompressibleSlipWallBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleSlipWallBoundaryField
 * \brief Inviscid slip-wall boundary condition for compressible DG states.
 *
 * \details
 * This boundary condition keeps density and total energy unchanged and
 * reflects only the wall-normal component of momentum:
 *
 * \f[
 *     \rho^+ = \rho^-,
 *     \qquad
 *     E^+ = E^-,
 * \f]
 *
 * \f[
 *     \rho\mathbf{U}^+
 *     =
 *     \rho\mathbf{U}^-
 *     -
 *     2\left(\mathbf{n}\cdot\rho\mathbf{U}^-\right)\mathbf{n}.
 * \f]
 *
 * Numerically this is the same reflection used by the symmetry boundary
 * condition, but it is intended for physical wall patches of type \c wall.
 */
class dgCompressibleSlipWallBoundaryField
:
    public dgCompressibleBoundaryField
{
public:

    TypeName("compressible::slipWall");

    dgCompressibleSlipWallBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    virtual ~dgCompressibleSlipWallBoundaryField() = default;

    virtual void updateGhostState
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        scalar& rhoPlus,
        vector& rhoUPlus,
        scalar& EPlus
    ) const override;

    virtual void updateBCValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        scalar& rhoBC,
        vector& rhoUBC,
        scalar& EBC
    ) const override;

    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
