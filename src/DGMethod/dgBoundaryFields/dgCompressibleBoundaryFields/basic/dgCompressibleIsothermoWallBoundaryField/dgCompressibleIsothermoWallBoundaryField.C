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
    velocityWeakEnforcement_(false)
{
    TValue_ = dict.get<scalar>("TValue");
    rampTime_ = dict.lookupOrDefault<scalar>("rampTime", 0);
    velocityWeakEnforcement_ =
        dict.lookupOrDefault<bool>("velocityWeakEnforcement", false);

    if (rampTime_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'rampTime' must be non-negative, but got "
            << rampTime_ << nl
            << exit(FatalIOError);
    }
}


scalar dgCompressibleIsothermoWallBoundaryField::rampFraction() const
{
    if (rampTime_ <= SMALL)
    {
        return scalar(1);
    }

    const scalar timeValue = dgMesh_.getFvMesh().time().value();

    return min(max(timeValue/rampTime_, scalar(0)), scalar(1));
}


void dgCompressibleIsothermoWallBoundaryField::updateGhostState
(
    const label,
    const label,
    const label,
    const vector&,
    const scalar rhoMinus,
    const vector& rhoUMinus,
    const scalar EMinus,
    scalar& rhoPlus,
    vector& rhoUPlus,
    scalar& EPlus
) const
{
    const scalar alpha = rampFraction();
    rhoPlus = rhoMinus;
    rhoUPlus =
        velocityWeakEnforcement_
      ? (scalar(1) - alpha)*rhoUMinus
      : (scalar(1) - scalar(2)*alpha)*rhoUMinus;
    EPlus = EMinus;
}


void dgCompressibleIsothermoWallBoundaryField::updateBCValue
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
    const vector UBC = (scalar(1) - alpha)*UMinus;
    const scalar TBC = (scalar(1) - alpha)*TMinus + alpha*TValue_;
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


void dgCompressibleIsothermoWallBoundaryField::correctSelfDiffusionFlux
(
    const label,
    const label,
    const label,
    const vector& n,
    vector& massDiffFlux
) const
{
    if (thermo_.selfDiffusion())
    {
        massDiffFlux -= (massDiffFlux & n)*n;
    }
}


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
