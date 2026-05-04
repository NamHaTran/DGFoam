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

#include "dgCompressibleIsothermoWallBoundaryField.H"
#include "dgCompressibleTemperatureJumpModel.H"
#include "dgCompressibleVelocitySlipModel.H"
#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleIsothermoWallBoundaryField, 0);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressibleIsothermoWallBoundaryField,
    dictionary
);

/**
 * \brief Construct an isothermal compressible wall boundary condition.
 *
 * \details
 * The constructor reads the prescribed wall temperature, optional ramping and
 * velocity-enforcement controls, and constructs any requested run-time
 * selectable velocity-slip or temperature-jump model from its coefficient
 * sub-dictionary.
 */
dgCompressibleIsothermoWallBoundaryField::
dgCompressibleIsothermoWallBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict),
    TValue_(Zero),
    rampTime_(Zero),
    velocityWeakEnforcement_(false),
    enableUSlip_(false),
    enableTJump_(false),
    Uwall_(Zero),
    USlipModel_(nullptr),
    TJumpModel_(nullptr)
{
    TValue_ = dict.get<scalar>("TValue");
    rampTime_ = dict.lookupOrDefault<scalar>("rampTime", 0);
    velocityWeakEnforcement_ =
        dict.lookupOrDefault<bool>("velocityWeakEnforcement", false);
    enableUSlip_ = dict.lookupOrDefault<bool>("enableUSlip", false);
    enableTJump_ = dict.lookupOrDefault<bool>("enableTJump", false);
    Uwall_ = dict.lookupOrDefault<vector>("Uwall", vector::zero);

    if (enableUSlip_)
    {
        const dictionary& coeffDict =
            dict.found("USlipCoefficients")
          ? dict.subDict("USlipCoefficients")
          : dictionary::null;

        USlipModel_ = dgCompressibleVelocitySlipModel::New
        (
            dict.get<word>("USlipModel"),
            patch,
            dgMesh,
            thermo,
            coeffDict
        );
    }

    if (enableTJump_)
    {
        const dictionary& coeffDict =
            dict.found("TJumpCoefficients")
          ? dict.subDict("TJumpCoefficients")
          : dictionary::null;

        TJumpModel_ = dgCompressibleTemperatureJumpModel::New
        (
            dict.get<word>("TJumpModel"),
            patch,
            dgMesh,
            thermo,
            coeffDict
        );
    }

    if (rampTime_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'rampTime' must be non-negative, but got "
            << rampTime_ << nl
            << exit(FatalIOError);
    }
}


/**
 * \brief Compute the current interpolation factor for ramped wall data.
 *
 * \details
 * A non-positive ramp duration means the boundary target is applied
 * immediately.  Otherwise the returned value increases linearly with time and
 * is clamped to the interval \f$[0,1]\f$.
 */
scalar dgCompressibleIsothermoWallBoundaryField::rampFraction() const
{
    if (rampTime_ <= SMALL)
    {
        return scalar(1);
    }

    const scalar timeValue = dgMesh_.getFvMesh().time().value();

    return min(max(timeValue/rampTime_, scalar(0)), scalar(1));
}


/**
 * \brief Map an owner-cell local face to the OpenFOAM patch distance factor.
 *
 * \details
 * Slip and jump models use this coefficient for Robin-type interpolation
 * lengths.  The DG global face is converted to the patch-local index before
 * reading \c fvPatch::deltaCoeffs().
 */
scalar dgCompressibleIsothermoWallBoundaryField::patchDeltaCoeff
(
    const label cellID,
    const label faceLocalID
) const
{
    const label globalFaceID = dgMesh_.cells()[cellID]->faces()[faceLocalID];
    const label patchFaceID =
        dgMesh_.getLocalFaceID(globalFaceID, patch_.index());

    return patch_.deltaCoeffs()[patchFaceID];
}


/**
 * \brief Assemble reusable wall-model input at one boundary Gauss point.
 *
 * \details
 * The state stores owner conservative data, optional ghost placeholders,
 * geometry, normal distance coefficient, and owner-side velocity and
 * temperature gradients.  These gradients are passed through for Maxwell slip
 * and Smoluchowski temperature-jump corrections.
 */
dgCompressibleBoundaryState
dgCompressibleIsothermoWallBoundaryField::makeBoundaryState
(
    const label cellID,
    const label faceLocalID,
    const label localGauss,
    const vector& n,
    const scalar rhoMinus,
    const vector& rhoUMinus,
    const scalar EMinus,
    const scalar rhoPlus,
    const vector& rhoUPlus,
    const scalar EPlus
) const
{
    dgCompressibleBoundaryState state;

    state.cellID = cellID;
    state.faceLocalID = faceLocalID;
    state.localGauss = localGauss;
    state.n = n;
    state.normalDeltaCoeff = patchDeltaCoeff(cellID, faceLocalID);

    state.rhoMinus = rhoMinus;
    state.rhoUMinus = rhoUMinus;
    state.EMinus = EMinus;

    state.rhoPlus = rhoPlus;
    state.rhoUPlus = rhoUPlus;
    state.EPlus = EPlus;

    state.gradUMinus =
        thermo_.gradU().gaussFields()[cellID]
            .faceField().minusValueOnFace(faceLocalID, localGauss);
    state.gradTMinus =
        thermo_.gradT().gaussFields()[cellID]
            .faceField().minusValueOnFace(faceLocalID, localGauss);

    return state;
}


/**
 * \brief Return the target boundary velocity.
 *
 * \details
 * If slip is disabled, the prescribed \c Uwall_ value is used.  If slip is
 * enabled, \c Uwall_ is passed as the reference wall velocity to the selected
 * model, which returns the slip-corrected tangential boundary velocity.
 */
vector dgCompressibleIsothermoWallBoundaryField::wallVelocity
(
    const dgCompressibleBoundaryState& state,
    const vector& UMinus
) const
{
    if (enableUSlip_)
    {
        return USlipModel_->USlip(state, Uwall_);
    }

    return Uwall_;
}


/**
 * \brief Return the target boundary temperature.
 *
 * \details
 * If temperature jump is disabled, the prescribed \c TValue_ is imposed
 * directly.  If enabled, the selected jump model computes a wall-adjacent
 * temperature using \c TValue_ as the reference wall temperature.
 */
scalar dgCompressibleIsothermoWallBoundaryField::wallTemperature
(
    const dgCompressibleBoundaryState& state,
    const scalar
) const
{
    if (enableTJump_)
    {
        return TJumpModel_->TJump(state, TValue_);
    }

    return TValue_;
}


/**
 * \brief Compute the conservative ghost state for the isothermal wall.
 *
 * \details
 * The owner conservative state is converted to primitive form, optional
 * slip/jump models provide target wall velocity and temperature, and ramping
 * blends those targets from the owner state.  Strong velocity enforcement uses
 * \f$\mathbf{U}^+ = 2\mathbf{U}_b - \mathbf{U}^-\f$; weak enforcement uses
 * \f$\mathbf{U}^+ = \mathbf{U}_b\f$ directly.
 */
void dgCompressibleIsothermoWallBoundaryField::updateConservativeGhostState
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
) const
{
    const scalar rhoMinusSafe = max(rhoMinus, SMALL);
    const vector UMinus = rhoUMinus/rhoMinusSafe;
    const scalar kMinus = 0.5*magSqr(UMinus);
    const scalar heMinus = EMinus/rhoMinusSafe - kMinus;
    const scalar TMinus =
        thermo_.calcTemperatureFromRhoHe(cellID, rhoMinusSafe, heMinus);
    const scalar alpha = rampFraction();

    const dgCompressibleBoundaryState state =
        makeBoundaryState
        (
            cellID,
            faceLocalID,
            localGauss,
            n,
            rhoMinus,
            rhoUMinus,
            EMinus
        );

    const vector UTarget = wallVelocity(state, UMinus);
    const scalar TTarget = wallTemperature(state, TMinus);

    const vector UBC = (scalar(1) - alpha)*UMinus + alpha*UTarget;
    const scalar TBC = (scalar(1) - alpha)*TMinus + alpha*TTarget;

    rhoPlus = rhoMinus;

    const vector UPlus =
        velocityWeakEnforcement_
      ? UBC
      : scalar(2)*UBC - UMinus;

    const scalar hePlus = thermo_.calcHeFromRhoT(rhoMinusSafe, TBC);
    rhoUPlus = rhoMinusSafe*UPlus;

    if (thermo_.heIsInternalEnergy())
    {
        EPlus = rhoMinusSafe*(hePlus + 0.5*magSqr(UPlus));
    }
    else
    {
        const scalar pPlus = thermo_.eos().calcPFromRhoT(rhoMinusSafe, TBC);
        EPlus = rhoMinusSafe*(hePlus + 0.5*magSqr(UPlus)) - pPlus;
    }
}


/**
 * \brief Compute the conservative value stored as the boundary state.
 *
 * \details
 * This method stores the ramped boundary state itself: density is copied from
 * the owner, momentum uses \f$\rho^-\mathbf{U}_b\f$, and energy is reconstructed
 * from the thermodynamic model using the ramped boundary temperature.
 */
void dgCompressibleIsothermoWallBoundaryField::updatePrimitiveBCValue
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
) const
{
    const scalar rhoMinusSafe = max(rhoMinus, SMALL);
    const vector UMinus = rhoUMinus/rhoMinusSafe;
    const scalar kMinus = 0.5*magSqr(UMinus);
    const scalar heMinus = EMinus/rhoMinusSafe - kMinus;
    const scalar TMinus =
        thermo_.calcTemperatureFromRhoHe(cellID, rhoMinusSafe, heMinus);
    const scalar alpha = rampFraction();

    const dgCompressibleBoundaryState state =
        makeBoundaryState
        (
            cellID,
            faceLocalID,
            localGauss,
            n,
            rhoMinus,
            rhoUMinus,
            EMinus
        );

    const vector UTarget = wallVelocity(state, UMinus);
    const scalar TTarget = wallTemperature(state, TMinus);

    const vector UBC = (scalar(1) - alpha)*UMinus + alpha*UTarget;
    const scalar TBC = (scalar(1) - alpha)*TMinus + alpha*TTarget;
    const scalar rhoBCSafe = max(rhoMinus, SMALL);
    const scalar heBC = thermo_.calcHeFromRhoT(rhoBCSafe, TBC);

    rhoBC = rhoMinus;
    rhoUBC = rhoBCSafe*UBC;

    if (thermo_.heIsInternalEnergy())
    {
        EBC = rhoBCSafe*(heBC + 0.5*magSqr(UBC));
    }
    else
    {
        const scalar pBC = thermo_.eos().calcPFromRhoT(rhoBCSafe, TBC);
        EBC = rhoBCSafe*(heBC + 0.5*magSqr(UBC)) - pBC;
    }
}


/**
 * \brief Apply wall flux correction for self-diffusion mode.
 *
 * \details
 * When self-diffusion is enabled, the full diffusive flux through the wall is
 * set to zero.
 */
void dgCompressibleIsothermoWallBoundaryField::correctFlux
(
    const label,
    const label,
    const label,
    const vector& n,
    vector& flux
) const
{
    if (thermo_.selfDiffusion())
    {
        flux -= (flux & n)*n;
    }
}


/**
 * \brief Validate that the condition is attached to an OpenFOAM wall patch.
 */
void dgCompressibleIsothermoWallBoundaryField::checkPatchType() const
{
    if (this->patch_.type() != "wall")
    {
        FatalErrorInFunction
            << "Boundary condition " << this->type() << " can only "
            << "be applied to patch type:\n"
            << "    wall\n"
            << "but patch " << this->patch_.name()
            << " is of type " << this->patch_.type() << nl
            << exit(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
