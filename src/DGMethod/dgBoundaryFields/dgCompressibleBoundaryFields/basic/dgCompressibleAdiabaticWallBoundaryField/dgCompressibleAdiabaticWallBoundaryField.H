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

#ifndef Foam_dgCompressibleAdiabaticWallBoundaryField_H
#define Foam_dgCompressibleAdiabaticWallBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class dgCompressibleAdiabaticWallBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleAdiabaticWallBoundaryField
 * \brief Adiabatic no-slip wall condition for compressible DG fields.
 *
 * \details
 * Reconstructs the interior primitive state and imposes a no-slip wall with
 * zero normal gradients for pressure and temperature:
 *
 * \f[
 *     \mathbf{U}^+ = \mathbf{0}, \qquad
 *     T^+ = T^-, \qquad
 *     p^+ = p^- .
 * \f]
 *
 * The ghost conservative state is then recovered from
 * \f$(p^+, T^+, \mathbf{U}^+)\f$ via the active thermodynamic model.
 */
class dgCompressibleAdiabaticWallBoundaryField
:
    public dgCompressibleBoundaryField
{
public:

    TypeName("compressible::adiabaticWall");

    dgCompressibleAdiabaticWallBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressibleAdiabaticWallBoundaryField() = default;

    /**
     * \brief Build the adiabatic no-slip ghost state.
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
    virtual void updateConservativeGhostState
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

    virtual void updatePrimitiveBCValue
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

    virtual void correctGradient
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        vector& grad
    ) const override;

    virtual void correctFlux
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        vector& flux
    ) const override;

    /** \brief Check that the boundary condition is applied on a wall patch. */
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
