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

#include "dgCompressibleMaxwellVelocitySlip.H"
#include "addToRunTimeSelectionTable.H"
#include "dgThermoConservative.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleMaxwellVelocitySlip, 0);
addToRunTimeSelectionTable
(
    dgCompressibleVelocitySlipModel,
    dgCompressibleMaxwellVelocitySlip,
    dictionary
);

dgCompressibleMaxwellVelocitySlip::dgCompressibleMaxwellVelocitySlip
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleVelocitySlipModel(patch, dgMesh, thermo, dict),
    accommodationCoeff_(dict.lookupOrDefault<scalar>("accommodationCoeff", 1)),
    thermalCreep_(dict.lookupOrDefault<Switch>("thermalCreep", true))
{
    if (accommodationCoeff_ <= SMALL || accommodationCoeff_ > 2)
    {
        FatalIOErrorInFunction(dict)
            << "Entry accommodationCoeff must satisfy "
            << "0 < accommodationCoeff <= 2"
            << exit(FatalIOError);
    }
}


vector dgCompressibleMaxwellVelocitySlip::USlip
(
    const dgCompressibleBoundaryState& state,
    const vector& Uwall
) const
{
    /*
     * Maxwell velocity slip:
     *
     *   C1 = sqrt((rho/p) pi/2) (2 - alphaU)/alphaU
     *   L_U = C1 nu
     *
     * The normal velocity-gradient contribution is treated implicitly with
     * the same sign convention as the Smoluchowski temperature jump:
     *
     *   grad_n(U_t) ~= deltaCoeff (U_t,b - U_t,minus)
     *   U_t,b = U_ref,t - L_U deltaCoeff (U_t,b - U_t,minus)
     *
     * where
     *
     *   U_ref = U_wall - 3 nu/(4 T) P_t grad(T)
     *
     * and P_t = I - n n is the tangential projector.  Solving the scalar
     * Robin equation in the tangential plane gives
     *
     *   U_t,b = (U_ref,t + deltaCoeff L_U U_minus,t)
     *           /(1 + deltaCoeff L_U),
     *
     * with the normal component fixed to the wall-normal velocity.  No
     * explicit grad(U) is required.
     */
    const scalar rhoMinus = max(state.rhoMinus, SMALL);
    const vector UMinus = state.rhoUMinus/rhoMinus;
    const scalar heMinus =
        state.EMinus/rhoMinus - scalar(0.5)*magSqr(UMinus);
    const scalar TMinus =
        thermo_.calcTemperatureFromRhoHe(state.cellID, rhoMinus, heMinus);
    const scalar pMinus =
        thermo_.calcPressureFromRhoHe(state.cellID, rhoMinus, heMinus);

    const scalar muMinus = thermo_.transport().calcMu(TMinus);
    const scalar psiMinus = rhoMinus/max(pMinus, SMALL);
    const scalar nuMinus = muMinus/rhoMinus;
    const scalar C1 =
        sqrt(psiMinus*constant::mathematical::piByTwo)
       *(scalar(2) - accommodationCoeff_)/accommodationCoeff_;

    const tensor tangentialProjector = tensor::I - state.n*state.n;

    vector refValue = Uwall;

    if (thermalCreep_)
    {
        refValue -=
            scalar(3)*nuMinus/(scalar(4)*max(TMinus, SMALL))
           *transform(tangentialProjector, state.gradTMinus);
    }

    const scalar slipLength = C1*nuMinus;
    const scalar valueFraction =
        scalar(1)/(scalar(1) + state.normalDeltaCoeff*slipLength);

    const vector UwallNormal = (Uwall & state.n)*state.n;
    const vector refTangential = transform(tangentialProjector, refValue);
    const vector UMinusTangential = transform(tangentialProjector, UMinus);

    return
        UwallNormal
      + valueFraction*refTangential
      + (scalar(1) - valueFraction)*UMinusTangential;
}

} // End namespace Foam

// ************************************************************************* //
