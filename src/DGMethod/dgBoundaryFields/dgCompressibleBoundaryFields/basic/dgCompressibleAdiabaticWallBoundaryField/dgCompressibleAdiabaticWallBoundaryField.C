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

#include "dgCompressibleAdiabaticWallBoundaryField.H"
#include "dgCompressibleVelocitySlipModel.H"
#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleAdiabaticWallBoundaryField, 0);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressibleAdiabaticWallBoundaryField,
    dictionary
);

/**
 * \brief Construct an adiabatic compressible wall boundary condition.
 *
 * \details
 * The constructor keeps the thermal behavior fixed to adiabatic and only reads
 * optional velocity-wall data.  If \c enableUSlip is true, the named
 * run-time-selectable velocity-slip model is constructed from the optional
 * \c USlipCoefficients sub-dictionary.
 */
dgCompressibleAdiabaticWallBoundaryField::
dgCompressibleAdiabaticWallBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict),
    enableUSlip_(dict.lookupOrDefault<bool>("enableUSlip", false)),
    Uwall_(dict.lookupOrDefault<vector>("Uwall", vector::zero)),
    USlipModel_(nullptr)
{
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
}


/**
 * \brief Map an owner-cell local face to the OpenFOAM patch distance factor.
 *
 * \details
 * Slip models use this coefficient for Robin-type blending between owner and
 * wall states.  The DG face numbering is first converted to a patch-local face
 * index before accessing \c fvPatch::deltaCoeffs().
 */
scalar dgCompressibleAdiabaticWallBoundaryField::patchDeltaCoeff
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
 * The adiabatic wall provides the velocity gradient for Maxwell slip
 * curvature/stress corrections, but sets \c gradTMinus to zero.  This keeps
 * temperature-jump physics and thermal-creep forcing out of the adiabatic
 * boundary while still allowing velocity slip.
 */
dgCompressibleBoundaryState
dgCompressibleAdiabaticWallBoundaryField::makeBoundaryState
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
    state.gradTMinus = Zero;

    return state;
}


/**
 * \brief Return the wall velocity imposed by the adiabatic boundary.
 *
 * \details
 * Without a slip model this is the prescribed \c Uwall_ value.  With a slip
 * model, \c Uwall_ is used as the wall reference velocity and the model returns
 * the tangentially corrected boundary velocity.
 */
vector dgCompressibleAdiabaticWallBoundaryField::wallVelocity
(
    const dgCompressibleBoundaryState& state,
    const vector&
) const
{
    if (enableUSlip_)
    {
        return USlipModel_->USlip(state, Uwall_);
    }

    return Uwall_;
}


/**
 * \brief Compute the conservative ghost state for the adiabatic wall.
 *
 * \details
 * The owner conservative variables are converted to \f$\mathbf{U}^-\f$,
 * \f$T^-\f$, and \f$p^-\f$.  The ghost state then uses \f$p^-=p^+\f$ and
 * \f$T^-=T^+\f$, with only the velocity replaced by the wall/slip velocity
 * before converting back to conservative form.
 */
void dgCompressibleAdiabaticWallBoundaryField::updateConservativeGhostState
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
    const scalar pMinus =
        thermo_.calcPressureFromRhoHe(cellID, rhoMinusSafe, heMinus);

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

    primitiveToConservative
    (
        pMinus,
        TMinus,
        wallVelocity(state, UMinus),
        rhoPlus,
        rhoUPlus,
        EPlus
    );
}


/**
 * \brief Compute the conservative value stored as the boundary state.
 *
 * \details
 * This mirrors \c updateConservativeGhostState: pressure and temperature are
 * reconstructed from the owner state, temperature is not jumped, and velocity
 * is supplied by \c wallVelocity().
 */
void dgCompressibleAdiabaticWallBoundaryField::updatePrimitiveBCValue
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
    const scalar pMinus =
        thermo_.calcPressureFromRhoHe(cellID, rhoMinusSafe, heMinus);

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

    primitiveToConservative
    (
        pMinus,
        TMinus,
        wallVelocity(state, UMinus),
        rhoBC,
        rhoUBC,
        EBC
    );
}


/**
 * \brief Remove the normal component from a boundary gradient.
 *
 * \details
 * This implements the adiabatic zero-normal-gradient correction by projecting
 * the supplied gradient into the tangent plane of the wall.
 */
void dgCompressibleAdiabaticWallBoundaryField::correctGradient
(
    const label,
    const label,
    const label,
    const vector& n,
    vector& grad
) const
{
    grad -= (grad & n)*n;
}


/**
 * \brief Apply wall flux correction for self-diffusion mode.
 *
 * \details
 * When self-diffusion is enabled, the full diffusive flux through the wall is
 * set to zero, enforcing the impermeable adiabatic boundary.
 */
void dgCompressibleAdiabaticWallBoundaryField::correctFlux
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
void dgCompressibleAdiabaticWallBoundaryField::checkPatchType() const
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
