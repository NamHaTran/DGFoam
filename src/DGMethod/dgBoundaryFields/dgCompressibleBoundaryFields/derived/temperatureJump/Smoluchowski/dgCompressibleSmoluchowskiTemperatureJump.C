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

#include "dgCompressibleSmoluchowskiTemperatureJump.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "dgThermoConservative.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleSmoluchowskiTemperatureJump, 0);
addToRunTimeSelectionTable
(
    dgCompressibleTemperatureJumpModel,
    dgCompressibleSmoluchowskiTemperatureJump,
    dictionary
);

dgCompressibleSmoluchowskiTemperatureJump::
dgCompressibleSmoluchowskiTemperatureJump
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleTemperatureJumpModel(patch, dgMesh, thermo, dict),
    accommodationCoeff_(dict.lookupOrDefault<scalar>("accommodationCoeff", 1)),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.4)),
    useThermoGamma_(dict.lookupOrDefault<bool>("useThermoGamma", true))
{
    if (accommodationCoeff_ <= SMALL || accommodationCoeff_ > 2)
    {
        FatalIOErrorInFunction(dict)
            << "Entry accommodationCoeff must satisfy "
            << "0 < accommodationCoeff <= 2"
            << exit(FatalIOError);
    }
}


scalar dgCompressibleSmoluchowskiTemperatureJump::TJump
(
    const dgCompressibleBoundaryState& state,
    const scalar Twall
) const
{
    /*
     * Smoluchowski temperature jump:
     *
     *   L_T = nu sqrt((rho/p) pi/2)
     *       * 2 gamma/(Pr (gamma + 1))
     *       * (2 - alphaT)/alphaT
     *
     * with alphaT = accommodationCoeff_.  The Robin interpolation below is
     * equivalent to T_b = (T_wall + deltaCoeff L_T T_minus)
     *                    /(1 + deltaCoeff L_T).
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
    const scalar Pr = thermo_.transport().calcPr(TMinus);
    const scalar Cp = thermo_.thermo().calcCp(TMinus);
    const scalar Cv = thermo_.thermo().calcCv(TMinus);
    const scalar gamma =
        useThermoGamma_ ? thermo_.thermo().calcGamma(Cp, Cv) : gamma_;
    const scalar psiMinus = rhoMinus/max(pMinus, SMALL);
    const scalar nuMinus = muMinus/rhoMinus;

    const scalar jumpLength =
        nuMinus
       *sqrt(psiMinus*constant::mathematical::piByTwo)
       *scalar(2)*gamma/(Pr*(gamma + scalar(1)))
       *(scalar(2) - accommodationCoeff_)/accommodationCoeff_;

    const scalar valueFraction =
        scalar(1)/(scalar(1) + state.normalDeltaCoeff*jumpLength);

    return valueFraction*Twall + (scalar(1) - valueFraction)*TMinus;
}

} // End namespace Foam

// ************************************************************************* //
