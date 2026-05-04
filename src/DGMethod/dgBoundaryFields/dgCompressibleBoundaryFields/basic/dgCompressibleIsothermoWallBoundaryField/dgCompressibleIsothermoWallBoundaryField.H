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

#include "autoPtr.H"
#include "dgCompressibleBoundaryField.H"
#include "dgCompressibleBoundaryState.H"

namespace Foam
{

class dgCompressibleTemperatureJumpModel;
class dgCompressibleVelocitySlipModel;

/*---------------------------------------------------------------------------*\
           Class dgCompressibleIsothermoWallBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleIsothermoWallBoundaryField
 * \brief Isothermal no-slip wall condition for compressible DG fields.
 *
 * \details
 * This wall model reconstructs the interior primitive pressure and then
 * imposes a wall state with prescribed temperature and wall velocity.  The
 * optional \c rampTime entry linearly ramps the target wall values from the
 * local owner state:
 *
 * \f[
 *     \mathbf{U}_b = \mathbf{U}_w, \qquad
 *     T^+ = T_w, \qquad
 *     p^+ = p^- .
 * \f]
 *
 * When \c enableUSlip is true, \f$\mathbf{U}_b\f$ is replaced by the selected
 * velocity-slip model value.  When \c enableTJump is true, \f$T_w\f$ is
 * replaced by the selected temperature-jump model value.  Unlike the
 * adiabatic wall, the owner temperature gradient is passed through the
 * reusable boundary state, so Maxwell thermal creep can use it.
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
 *
 * Required/optional dictionary entries:
 * \verbatim
 * TValue                    300;       // Required fixed wall temperature
 * rampTime                  0;         // Optional linear ramp duration
 * velocityWeakEnforcement   false;     // Optional use UBC directly as ghost U
 * enableUSlip               false;     // Optional velocity-slip correction
 * USlipModel                Maxwell;   // Required only when enableUSlip true
 * USlipCoefficients         { ... }    // Optional model coefficients
 * enableTJump               false;     // Optional temperature-jump correction
 * TJumpModel                Smoluchowski;
 * TJumpCoefficients         { ... }    // Optional model coefficients
 * Uwall                     (0 0 0);   // Optional prescribed wall velocity
 * \endverbatim
 */
class dgCompressibleIsothermoWallBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    /** \brief Prescribed wall temperature. */
    scalar TValue_;

    /**
     * \brief Linear ramp duration for fixed wall values.
     *
     * \details
     * Values less than or equal to \c SMALL disable ramping and impose the
     * target wall state immediately.
     */
    scalar rampTime_;

    /**
     * \brief Select weak or mirrored velocity enforcement for the ghost state.
     *
     * \details
     * When true, the ghost velocity is the target boundary velocity.  When
     * false, the ghost velocity is mirrored about the boundary value as
     * \f$\mathbf{U}^+ = 2\mathbf{U}_b - \mathbf{U}^-\f$.
     */
    bool velocityWeakEnforcement_;

    /** \brief Enable run-time selectable velocity-slip model. */
    bool enableUSlip_;

    /** \brief Enable run-time selectable temperature-jump model. */
    bool enableTJump_;

    /** \brief Prescribed wall velocity used as no-slip/slip reference. */
    vector Uwall_;

    /** \brief Optional velocity-slip model, allocated only when enabled. */
    autoPtr<dgCompressibleVelocitySlipModel> USlipModel_;

    /** \brief Optional temperature-jump model, allocated only when enabled. */
    autoPtr<dgCompressibleTemperatureJumpModel> TJumpModel_;

    /**
     * \brief Return the current linear ramp fraction.
     *
     * \return Clamped value in the interval \f$[0,1]\f$.
     */
    scalar rampFraction() const;

    /**
     * \brief Return the finite-volume normal distance coefficient for a face.
     *
     * \param cellID Owner cell index.
     * \param faceLocalID Local face index in the owner DG cell.
     * \return Patch \c deltaCoeffs value for the corresponding boundary face.
     */
    scalar patchDeltaCoeff
    (
        const label cellID,
        const label faceLocalID
    ) const;

    /**
     * \brief Gather all quantities required by reusable slip/jump models.
     *
     * \details
     * The state includes owner conservative variables, optional ghost
     * placeholders, wall normal, patch distance coefficient, and owner-side
     * velocity/temperature gradients at the same Gauss point.
     */
    dgCompressibleBoundaryState makeBoundaryState
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        const scalar rhoPlus = Zero,
        const vector& rhoUPlus = vector::zero,
        const scalar EPlus = Zero
    ) const;

    /**
     * \brief Compute the target wall velocity.
     *
     * \details
     * Returns \c Uwall_ when slip is disabled, otherwise delegates to the
     * selected velocity-slip model.
     */
    vector wallVelocity
    (
        const dgCompressibleBoundaryState& state,
        const vector& UMinus
    ) const;

    /**
     * \brief Compute the target wall temperature.
     *
     * \details
     * Returns \c TValue_ when temperature jump is disabled, otherwise delegates
     * to the selected temperature-jump model using \c TValue_ as the wall
     * reference temperature.
     */
    scalar wallTemperature
    (
        const dgCompressibleBoundaryState& state,
        const scalar TMinus
    ) const;

public:

    TypeName("compressible::isothermoWall");

    /**
     * \brief Construct from patch, DG mesh, thermodynamics, and dictionary.
     *
     * \details
     * Reads fixed wall temperature, optional ramp/velocity enforcement
     * settings, and constructs the optional slip/jump models.
     */
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
     * \brief Build the conservative ghost state for the isothermal wall.
     *
     * \details
     * Reconstructs \f$T^-\f$ from the owner state, evaluates target
     * velocity/temperature including optional slip/jump corrections, applies
     * the ramp fraction, and converts the ghost primitive state to
     * conservative form.
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

    /**
     * \brief Build the conservative value used as the primitive boundary state.
     *
     * \details
     * Stores the ramped target wall state directly, without the mirrored ghost
     * velocity used for strong enforcement.
     */
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

    /**
     * \brief Correct diffusive/self-diffusion fluxes at the wall.
     *
     * \details
     * When self diffusion is active, the boundary flux is forced to zero.
     */
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
