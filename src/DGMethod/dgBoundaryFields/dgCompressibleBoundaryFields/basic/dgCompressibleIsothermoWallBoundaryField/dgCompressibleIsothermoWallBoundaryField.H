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

#ifndef Foam_dgCompressibleIsothermoWallBoundaryField_H
#define Foam_dgCompressibleIsothermoWallBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class dgCompressibleIsothermoWallBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleIsothermoWallBoundaryField
 * \brief Isothermal no-slip wall condition for compressible DG fields.
 *
 * \details
 * This wall model reconstructs the interior primitive pressure and then
 * imposes a wall state with prescribed temperature and zero velocity:
 *
 * \f[
 *     \mathbf{U}^+ = \mathbf{0}, \qquad
 *     T^+ = T_w, \qquad
 *     p^+ = p^- .
 * \f]
 *
 * The interior pressure is recovered from the conservative interior state via
 *
 * \f[
 *     \mathbf{U}^- = \frac{\rho\mathbf{U}^-}{\rho^-}, \qquad
 *     k^- = \frac{1}{2}\left|\mathbf{U}^-\right|^2,
 * \f]
 *
 * \f[
 *     h_e^- = \frac{E^-}{\rho^-} - k^-,
 * \qquad
 *     T^- = T(h_e^-),
 * \qquad
 *     p^- = p(\rho^-, T^-).
 * \f]
 *
 * The ghost conservative state is then obtained from
 * \f$(p^+, T^+, \mathbf{U}^+)\f$ using the active thermodynamic model.
 */
class dgCompressibleIsothermoWallBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    /** \brief Prescribed wall temperature. */
    scalar TValue_;

public:

    TypeName("compressible::isothermoWall");

    dgCompressibleIsothermoWallBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressibleIsothermoWallBoundaryField() = default;

    /**
     * \brief Build the isothermal no-slip ghost state.
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

    /** \brief Check that the boundary condition is applied on a wall patch. */
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
