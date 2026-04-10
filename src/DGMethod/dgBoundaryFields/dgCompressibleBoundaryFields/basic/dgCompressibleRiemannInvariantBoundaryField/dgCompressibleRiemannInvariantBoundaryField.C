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

#include "dgCompressibleRiemannInvariantBoundaryField.H"
#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleRiemannInvariantBoundaryField, 0);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressibleRiemannInvariantBoundaryField,
    dictionary
);

namespace
{

template<class Type>
Type readInfOrValue(const dictionary& dict, const word& infName, const word& valueName)
{
    if (dict.found(infName))
    {
        return dict.get<Type>(infName);
    }

    return dict.get<Type>(valueName);
}

} // End anonymous namespace


dgCompressibleRiemannInvariantBoundaryField::
dgCompressibleRiemannInvariantBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict),
    pInf_(Zero),
    TInf_(Zero),
    UInf_(Zero),
    rhoInf_(Zero),
    gammaInf_(Zero),
    cInf_(Zero),
    sInf_(Zero),
    fallbackToZeroGradientIfTangential_(true),
    tangentialAngleTolDeg_(5.0),
    tangentialMinSpeed_(1e-6)
{
    const dictionary& coeffDict = dict.subDict("compressibleRiemannInvariantCoeff");

    pInf_ = readInfOrValue<scalar>(coeffDict, "pInf", "pValue");
    TInf_ = readInfOrValue<scalar>(coeffDict, "TInf", "TValue");
    UInf_ = readInfOrValue<vector>(coeffDict, "UInf", "UValue");

    rhoInf_ = thermo_.eos().calcRhoFromPT(pInf_, TInf_);
    const scalar heInf = thermo_.calcHeFromRhoT(rhoInf_, TInf_);

    const scalar CpInf = thermo_.thermo().calcCp(TInf_);
    const scalar CvInf = thermo_.thermo().calcCv(TInf_);
    gammaInf_ = thermo_.thermo().calcGamma(CpInf, CvInf);
    cInf_ = thermo_.calcSpeedOfSoundFromRhoHe(rhoInf_, heInf);

    fallbackToZeroGradientIfTangential_ =
        coeffDict.lookupOrDefault<bool>
        (
            "fallbackToZeroGradientIfTangential",
            true
        );

    tangentialAngleTolDeg_ =
        coeffDict.lookupOrDefault<scalar>
        (
            "tangentialAngleTolDeg",
            5.0
        );

    tangentialMinSpeed_ =
        coeffDict.lookupOrDefault<scalar>
        (
            "tangentialMinSpeed",
            1e-6
        );

    tangentialAngleTolDeg_ =
        min(max(tangentialAngleTolDeg_, scalar(0)), scalar(90));
    tangentialMinSpeed_ = max(tangentialMinSpeed_, scalar(0));

    // This characteristic form assumes a calorically perfect-gas style
    // relation with one reference gamma.
    sInf_ = pInf_/pow(rhoInf_, gammaInf_);
}


void dgCompressibleRiemannInvariantBoundaryField::updateGhostState
(
    const label,
    const label,
    const label,
    const vector& n,
    const scalar rhoMinus,
    const vector& rhoUMinus,
    const scalar EMinus,
    scalar& rhoPlus,
    vector& rhoUPlus,
    scalar& EPlus
) const
{
    const scalar gammaMinusOne = gammaInf_ - 1.0;
    const scalar gammaInv = 1.0/gammaInf_;
    const scalar gammaMinusOneInv = 1.0/gammaMinusOne;

    const vector UMinus = rhoUMinus/rhoMinus;
    const scalar VnMinus = n & UMinus;
    const scalar VnInf = n & UInf_;
    const scalar magUMinus = mag(UMinus);

    if
    (
        fallbackToZeroGradientIfTangential_
     && magUMinus > tangentialMinSpeed_
    )
    {
        const scalar cosTheta = min(mag(VnMinus)/magUMinus, scalar(1));
        const scalar tangentialCosTol =
            Foam::sin(degToRad(tangentialAngleTolDeg_));

        if (cosTheta <= tangentialCosTol)
        {
            rhoPlus = rhoMinus;
            rhoUPlus = rhoUMinus;
            EPlus = EMinus;
            return;
        }
    }

    const scalar kMinus = 0.5*magSqr(UMinus);
    const scalar heMinus = EMinus/rhoMinus - kMinus;
    const scalar pMinus = thermo_.calcPressureFromRhoHe(rhoMinus, heMinus);
    const scalar cMinusLocal = sqrt(gammaInf_*pMinus/rhoMinus);
    const scalar machMinus = mag(VnMinus)/max(cMinusLocal, SMALL);

    scalar rPlus = Zero;
    scalar rMinus = Zero;

    if (VnMinus <= 0.0)
    {
        if (machMinus < 1.0)
        {
            rPlus = VnMinus + 2.0*cMinusLocal*gammaMinusOneInv;
            rMinus = VnInf - 2.0*cInf_*gammaMinusOneInv;
        }
        else
        {
            rPlus = VnInf + 2.0*cInf_*gammaMinusOneInv;
            rMinus = VnInf - 2.0*cInf_*gammaMinusOneInv;
        }
    }
    else
    {
        if (machMinus < 1.0)
        {
            rPlus = VnMinus + 2.0*cMinusLocal*gammaMinusOneInv;
            rMinus = VnInf - 2.0*cInf_*gammaMinusOneInv;
        }
        else
        {
            rPlus = VnMinus + 2.0*cMinusLocal*gammaMinusOneInv;
            rMinus = VnMinus - 2.0*cMinusLocal*gammaMinusOneInv;
        }
    }

    const scalar VnBC = 0.5*(rPlus + rMinus);
    const scalar cBC = 0.25*gammaMinusOne*(rPlus - rMinus);

    vector UBC(Zero);
    scalar sBC = Zero;

    if (VnMinus <= 0.0)
    {
        UBC = UInf_ + (VnBC - VnInf)*n;
        sBC = sInf_;
    }
    else
    {
        UBC = UMinus + (VnBC - VnMinus)*n;
        sBC = pMinus/pow(rhoMinus, gammaInf_);
    }

    rhoPlus = pow((cBC*cBC)/(gammaInf_*sBC), gammaMinusOneInv);
    const scalar pBC = rhoPlus*cBC*cBC*gammaInv;
    rhoUPlus = rhoPlus*UBC;

    if (thermo_.heIsInternalEnergy())
    {
        EPlus = pBC*gammaMinusOneInv + 0.5*rhoPlus*magSqr(UBC);
    }
    else
    {
        EPlus = pBC*gammaMinusOneInv + 0.5*rhoPlus*magSqr(UBC) - pBC;
    }
}


void dgCompressibleRiemannInvariantBoundaryField::checkPatchType() const
{
    if (this->patch_.type() != "patch")
    {
        FatalErrorInFunction
            << "Boundary condition " << this->type() << " can only "
            << "be applied to patch type:\n"
            << "    patch\n"
            << "but patch " << this->patch_.name()
            << " is of type " << this->patch_.type() << nl
            << exit(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
