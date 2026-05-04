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

#ifndef Foam_dgCompressibleMaxwellVelocitySlip_H
#define Foam_dgCompressibleMaxwellVelocitySlip_H

#include "Switch.H"
#include "dgCompressibleVelocitySlipModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dgCompressibleMaxwellVelocitySlip Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleMaxwellVelocitySlip
 * \brief Maxwell velocity-slip model for rarefied wall boundaries.
 *
 * \details
 * The model corrects the tangential wall velocity with optional thermal creep,
 * then solves the Maxwell normal-gradient slip relation implicitly with a
 * forward-difference approximation for the tangential velocity gradient.
 */
class dgCompressibleMaxwellVelocitySlip
:
    public dgCompressibleVelocitySlipModel
{
private:

    /// Tangential momentum accommodation coefficient.
    scalar accommodationCoeff_;

    /// Include the thermal-creep contribution.
    Switch thermalCreep_;

public:

    /// Run-time type name.
    TypeName("Maxwell");

    /// Construct from the patch, mesh, thermo model and coefficient dictionary.
    dgCompressibleMaxwellVelocitySlip
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /**
     * \brief Compute the Maxwell wall velocity slip.
     *
     * \param state Interior/exterior state at one boundary Gauss point.
     * \param Uwall Prescribed wall velocity.
     * \return Effective boundary velocity after the slip correction.
     */
    virtual vector USlip
    (
        const dgCompressibleBoundaryState& state,
        const vector& Uwall
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
