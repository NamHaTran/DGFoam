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

#ifndef Foam_dgCompressiblePressureMachTemperatureBoundaryField_H
#define Foam_dgCompressiblePressureMachTemperatureBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class dgCompressiblePressureMachTemperatureBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressiblePressureMachTemperatureBoundaryField
 * \brief Compressible DG boundary condition prescribed by pressure, temperature, and Mach.
 *
 * \details
 * Builds a constant ghost state from user-specified
 * \f$(p_b,\,T_b,\,\mathbf{M}_b)\f$ data.
 *
 * The Mach vector is converted into velocity through the local speed of sound:
 *
 * \f[
 *     \mathbf{U}_b = a_b \mathbf{M}_b,
 * \f]
 *
 * where \f$a_b\f$ is the thermodynamic speed of sound evaluated from the
 * prescribed boundary state. The ghost conservative state is then
 *
 * \f[
 *     (\rho^+,\,\rho\mathbf{U}^+,\,E^+)
 *     =
 *     \mathcal{C}\!\left(p_b,\,T_b,\,\mathbf{U}_b\right).
 * \f]
 *
 * The resulting state is constant over the whole patch and does not depend on
 * the interior state.
 */
class dgCompressiblePressureMachTemperatureBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    /** \brief Prescribed boundary density in conservative form. */
    scalar rhoValue_;

    /** \brief Prescribed boundary momentum in conservative form. */
    vector rhoUValue_;

    /** \brief Prescribed boundary total energy in conservative form. */
    scalar EValue_;

public:

    TypeName("compressible::pressureMachTemperature");

    dgCompressiblePressureMachTemperatureBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressiblePressureMachTemperatureBoundaryField() = default;

    /**
     * \brief Assign the precomputed constant ghost state.
     *
     * \param cellID Owner cell index.
     * \param faceLocalID Local face index on the owner cell.
     * \param localGauss Local Gauss-point index on the face.
     * \param n Outward unit normal of the owner cell.
     * \param rhoMinus Interior density.
     * \param rhoUMinus Interior momentum.
     * \param EMinus Interior total energy.
     * \param rhoPlus Ghost density.
     * \param rhoUPlus Ghost momentum.
     * \param EPlus Ghost total energy.
     */
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

    /** \brief Check that the boundary condition is applied on a supported patch type. */
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
