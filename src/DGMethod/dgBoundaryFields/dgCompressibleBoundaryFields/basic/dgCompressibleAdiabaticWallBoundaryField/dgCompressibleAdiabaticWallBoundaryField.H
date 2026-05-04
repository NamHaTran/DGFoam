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

#include "autoPtr.H"
#include "dgCompressibleBoundaryField.H"
#include "dgCompressibleBoundaryState.H"

namespace Foam
{

class dgCompressibleVelocitySlipModel;

/*---------------------------------------------------------------------------*\
           Class dgCompressibleAdiabaticWallBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleAdiabaticWallBoundaryField
 * \brief Adiabatic wall condition for compressible DG fields.
 *
 * \details
 * Reconstructs the owner-side primitive state at each boundary Gauss point and
 * imposes an adiabatic wall.  The thermal condition is strictly zero heat
 * flux: the boundary temperature is the reconstructed owner temperature and no
 * temperature-jump model is applied by this boundary condition.
 *
 * \f[
 *     \mathbf{U}^+ = \mathbf{U}_w, \qquad
 *     T^+ = T^-, \qquad
 *     p^+ = p^- .
 * \f]
 *
 * The wall velocity defaults to zero.  When \c enableUSlip is true, the
 * selected \c USlipModel computes the tangential wall state.  Since the wall
 * remains adiabatic, the reusable boundary state passes \c gradTMinus as zero
 * to velocity-slip models, disabling thermal-creep contributions that depend
 * on the temperature gradient.
 *
 * Required/optional dictionary entries:
 * \verbatim
 * enableUSlip        false;       // Optional Maxwell-type slip correction
 * USlipModel         Maxwell;     // Required only when enableUSlip is true
 * USlipCoefficients  { ... }      // Optional model-specific coefficients
 * Uwall              (0 0 0);     // Optional prescribed wall velocity
 * \endverbatim
 *
 * The ghost conservative state is recovered from
 * \f$(p^-, T^-, \mathbf{U}_b)\f$ via the active thermodynamic model, where
 * \f$\mathbf{U}_b\f$ is either \c Uwall or the slip-model velocity.
 */
class dgCompressibleAdiabaticWallBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    /** \brief Enable the run-time selectable velocity-slip model. */
    bool enableUSlip_;

    /** \brief Prescribed wall velocity used as no-slip value or slip reference. */
    vector Uwall_;

    /** \brief Optional velocity-slip model, allocated only when enabled. */
    autoPtr<dgCompressibleVelocitySlipModel> USlipModel_;

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
     * \brief Gather all quantities required by reusable wall models.
     *
     * \details
     * The state stores conservative owner/ghost values, geometry, and the
     * owner-side velocity gradient.  For this adiabatic condition
     * \c gradTMinus is deliberately set to zero before passing the state to a
     * velocity-slip model.
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
     * \brief Compute the boundary velocity for the current Gauss point.
     *
     * \details
     * Returns \c Uwall_ when slip is disabled, otherwise delegates to
     * \c USlipModel_ using the preassembled boundary state.
     */
    vector wallVelocity
    (
        const dgCompressibleBoundaryState& state,
        const vector& UMinus
    ) const;

public:

    TypeName("compressible::adiabaticWall");

    /**
     * \brief Construct from patch, DG mesh, thermodynamics, and dictionary.
     *
     * \details
     * Reads \c enableUSlip, optional \c Uwall, and constructs the selected
     * velocity-slip model when requested.  Temperature-jump entries are ignored
     * by design for an adiabatic wall.
     */
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
     * \brief Build the conservative ghost state for the adiabatic wall.
     *
     * \details
     * Reconstructs \f$T^-\f$ and \f$p^-\f$ from the owner conservative state,
     * evaluates the wall/slip velocity, keeps \f$T_b=T^-\f$, and converts the
     * resulting primitive wall state to conservative ghost variables.
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
     * Uses the same target state as \c updateConservativeGhostState so that
     * strongly imposed and diagnostic boundary values remain consistent.
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
     * \brief Enforce zero normal gradient at an adiabatic wall.
     *
     * \details
     * Removes the normal component from the supplied gradient, leaving only the
     * tangential part.
     */
    virtual void correctGradient
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        vector& grad
    ) const override;

    /**
     * \brief Correct diffusive/self-diffusion fluxes at the wall.
     *
     * \details
     * When self diffusion is active, the boundary flux is forced to zero to
     * preserve the impermeable adiabatic-wall constraint.
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
