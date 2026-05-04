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

#ifndef Foam_dgCompressibleSmoluchowskiTemperatureJump_H
#define Foam_dgCompressibleSmoluchowskiTemperatureJump_H

#include "dgCompressibleTemperatureJumpModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class dgCompressibleSmoluchowskiTemperatureJump Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleSmoluchowskiTemperatureJump
 * \brief Smoluchowski temperature-jump model for rarefied wall boundaries.
 *
 * \details
 * The model blends the prescribed wall temperature and the interior
 * temperature using a Smoluchowski jump length based on viscosity,
 * pressure, Prandtl number, heat-capacity ratio and accommodation
 * coefficient.
 */
class dgCompressibleSmoluchowskiTemperatureJump
:
    public dgCompressibleTemperatureJumpModel
{
private:

    /// Thermal accommodation coefficient.
    scalar accommodationCoeff_;

    /// Constant heat-capacity ratio used when \c useThermoGamma_ is false.
    scalar gamma_;

    /// Use thermodynamic Cp/Cv for gamma when enabled.
    bool useThermoGamma_;

public:

    /// Run-time type name.
    TypeName("Smoluchowski");

    /// Construct from the patch, mesh, thermo model and coefficient dictionary.
    dgCompressibleSmoluchowskiTemperatureJump
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /**
     * \brief Compute the Smoluchowski wall temperature jump.
     *
     * \param state Interior/exterior state at one boundary Gauss point.
     * \param Twall Prescribed wall temperature.
     * \return Effective boundary temperature after the jump correction.
     */
    virtual scalar TJump
    (
        const dgCompressibleBoundaryState& state,
        const scalar Twall
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
