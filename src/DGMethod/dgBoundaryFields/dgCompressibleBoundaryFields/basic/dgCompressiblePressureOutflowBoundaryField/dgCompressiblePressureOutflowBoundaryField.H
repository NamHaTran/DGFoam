/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

#ifndef Foam_dgCompressiblePressureOutflowBoundaryField_H
#define Foam_dgCompressiblePressureOutflowBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class dgCompressiblePressureOutflowBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressiblePressureOutflowBoundaryField
 * \brief Static-pressure outflow condition for compressible DG fields.
 *
 * \details
 * The boundary inspects the interior normal Mach number
 *
 * \f[
 *     M_n^- = \frac{|\boldsymbol{n}\cdot\boldsymbol{U}^-|}{a^-},
 * \f]
 *
 * where
 *
 * \f[
 *     \boldsymbol{U}^- = \frac{\rho\boldsymbol{U}^-}{\rho^-},
 *     \qquad
 *     a^- = a(\rho^-, h_e^-),
 *     \qquad
 *     h_e^- = \frac{E^-}{\rho^-} - \frac{1}{2}\left|\boldsymbol{U}^-\right|^2.
 * \f]
 *
 * For subsonic outflow, density and momentum are extrapolated from the
 * interior while the pressure is replaced by a prescribed outlet value
 * \f$p_{\mathrm{out}}\f$. The boundary total energy is rebuilt from
 * \f$(p_{\mathrm{out}}, \rho^-, \boldsymbol{U}^-)\f$:
 *
 * \f[
 *     \rho^+ = \rho^-,
 *     \qquad
 *     \rho\boldsymbol{U}^+ = \rho\boldsymbol{U}^-,
 * \f]
 *
 * \f[
 *     T^+ = T(p_{\mathrm{out}}, \rho^-),
 *     \qquad
 *     h_e^+ = h_e(\rho^+, T^+),
 * \f]
 *
 * \f[
 *     E^+ =
 *     \begin{cases}
 *         \rho^+\!\left(h_e^+ + \frac{1}{2}\left|\boldsymbol{U}^-\right|^2\right),
 *         & \text{internal-energy formulation}, \\
 *         \rho^+\!\left(h_e^+ + \frac{1}{2}\left|\boldsymbol{U}^-\right|^2\right)
 *         - p_{\mathrm{out}},
 *         & \text{enthalpy formulation}.
 *     \end{cases}
 * \f]
 *
 * For supersonic outflow, all conservative variables are copied from the
 * interior state unchanged.
 */
class dgCompressiblePressureOutflowBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    /** \brief Prescribed static outlet pressure. */
    scalar pValue_;

    /** \brief Subsonic/supersonic split based on the interior normal Mach number. */
    scalar subsonicMachThreshold_;

public:

    TypeName("compressible::pressureOutflow");

    dgCompressiblePressureOutflowBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressiblePressureOutflowBoundaryField() = default;

    /**
     * \brief Build the outlet ghost state from the interior trace and target pressure.
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

    /** \brief Check that the boundary condition is applied on a standard patch. */
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
