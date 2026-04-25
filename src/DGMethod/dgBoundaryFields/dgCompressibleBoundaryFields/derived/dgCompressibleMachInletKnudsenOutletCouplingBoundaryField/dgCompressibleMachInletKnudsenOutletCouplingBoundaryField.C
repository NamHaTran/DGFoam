/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2026 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "dgCompressibleMachInletKnudsenOutletCouplingBoundaryField.H"
#include "dgCompressibleBoundaryManager.H"
#include "dgThermoConservative.H"
#include "transportLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug
(
    dgCompressibleMachInletKnudsenOutletCouplingBoundaryField,
    0
);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressibleMachInletKnudsenOutletCouplingBoundaryField,
    dictionary
);


dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict),
    position_(word::null),
    Kn_(Zero),
    MachIn_(Zero),
    Tin_(Zero),
    TOut_(Zero),
    pressureRatio_(1),
    refL_(1),
    pOut_(Zero),
    pIn_(Zero),
    direction_(Zero),
    usePatchNormalDirection_(true),
    coupledPatchPtr_(nullptr)
{
    const dictionary& coeffDict = dict.subDict
    (
        "compressibleMachInletKnudsenOutletCouplingCoeff"
    );

    position_ =
        dict.found("position")
      ? dict.get<word>("position")
      : coeffDict.get<word>("position");

    if (position_ != "inlet" && position_ != "outlet")
    {
        FatalIOErrorInFunction(dict)
            << "position must be either inlet or outlet, got "
            << position_ << exit(FatalIOError);
    }

    Tin_ = coeffDict.lookupOrDefault<scalar>("Tin", scalar(300));
    TOut_ = coeffDict.lookupOrDefault<scalar>("TOut", Tin_);
    MachIn_ = coeffDict.lookupOrDefault<scalar>("MachIn", scalar(0));
    pressureRatio_ =
        coeffDict.lookupOrDefault<scalar>("pressureRatio", scalar(1));
    Kn_ = coeffDict.lookupOrDefault<scalar>("Kn", scalar(0));
    refL_ = coeffDict.lookupOrDefault<scalar>("refL", scalar(1));

    if (coeffDict.found("direction"))
    {
        direction_ = coeffDict.get<vector>("direction");
        const scalar magDir = mag(direction_);

        if (magDir <= SMALL)
        {
            FatalIOErrorInFunction(dict)
                << "direction must be non-zero" << exit(FatalIOError);
        }

        direction_ /= magDir;
        usePatchNormalDirection_ = false;
    }

    if
    (
        position_ == "inlet"
     && coupling_
     && !coeffDict.found("pOutValue")
     && !coeffDict.found("Kn")
    )
    {
        pOut_ = Zero;
    }
    else
    {
        pOut_ = calcPout(coeffDict);
    }

    pIn_ = max(pressureRatio_*pOut_, scalar(SMALL));
}


scalar
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::calcPout
(
    const dictionary& coeffDict
) const
{
    if (coeffDict.found("pOutValue"))
    {
        return max(coeffDict.get<scalar>("pOutValue"), scalar(SMALL));
    }

    const scalar Kn = coeffDict.lookupOrDefault<scalar>("Kn", Kn_);
    const scalar L = coeffDict.lookupOrDefault<scalar>("refL", refL_);
    const scalar T =
        coeffDict.lookupOrDefault<scalar>("TOut", TOut_);

    if (Kn <= SMALL || L <= SMALL || T <= SMALL)
    {
        FatalIOErrorInFunction(dict_)
            << "Kn, refL, and TOut must be positive when "
            << "pOutValue is not supplied" << exit(FatalIOError);
    }

    const scalar mu = thermo_.transport().calcMu(T);

    return max
    (
        mu*sqrt(constant::mathematical::pi*thermo_.R()*T/scalar(2))/(Kn*L),
        scalar(SMALL)
    );
}


vector
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::inletVelocity
(
    const vector& n
) const
{
    const vector dir = usePatchNormalDirection_ ? -n : direction_;

    const scalar Cp = thermo_.thermo().calcCp(Tin_);
    const scalar Cv = thermo_.thermo().calcCv(Tin_);
    const scalar gamma = thermo_.thermo().calcGamma(Cp, Cv);
    const scalar a = thermo_.thermo().calcSpeedOfSound(Tin_, gamma);

    return MachIn_*a*dir;
}


void
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::updateInletState
(
    const vector& n,
    scalar& rhoPlus,
    vector& rhoUPlus,
    scalar& EPlus
) const
{
    const scalar pIn =
        coupledPatchPtr_ ? pressureRatio_*coupledPatchPtr_->pOut() : pIn_;

    primitiveToConservative
    (
        max(pIn, scalar(SMALL)),
        Tin_,
        inletVelocity(n),
        rhoPlus,
        rhoUPlus,
        EPlus
    );
}


scalar
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::pOut() const
{
    if (position_ != "outlet")
    {
        FatalErrorInFunction
            << "pOut() is only valid for outlet-position "
            << this->type() << " patches. Patch " << patch_.name()
            << " has position " << position_ << exit(FatalError);
    }

    return pOut_;
}


void dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::resolveCoupling
(
    const dgCompressibleBoundaryManager& manager
)
{
    dgCompressibleBoundaryField::resolveCoupling(manager);

    if (!coupling_)
    {
        return;
    }

    const dgCompressibleBoundaryField& coupled =
        manager.boundaryField(couplePatchName_);

    coupledPatchPtr_ = dynamic_cast
    <
        const dgCompressibleMachInletKnudsenOutletCouplingBoundaryField*
    >(&coupled);

    if (!coupledPatchPtr_)
    {
        FatalIOErrorInFunction(dict_)
            << "Patch " << patch_.name()
            << " couples to patch " << couplePatchName_
            << ", but that patch does not use boundary condition "
            << this->type() << exit(FatalIOError);
    }

    if (coupledPatchPtr_->position() == position_)
    {
        FatalIOErrorInFunction(dict_)
            << "Patch " << patch_.name()
            << " and couplePatch " << couplePatchName_
            << " both use position " << position_
            << ". Coupled patches must use opposite inlet/outlet positions."
            << exit(FatalIOError);
    }
}


void
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::updateConservativeGhostState
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
    if (position_ == "inlet")
    {
        updateInletState(n, rhoPlus, rhoUPlus, EPlus);
        return;
    }

    rhoPlus = rhoMinus;
    rhoUPlus = rhoUMinus;
    EPlus = EMinus;
}


void
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::updatePrimitiveBCValue
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
    updateConservativeGhostState
    (
        cellID,
        faceLocalID,
        localGauss,
        n,
        rhoMinus,
        rhoUMinus,
        EMinus,
        rhoBC,
        rhoUBC,
        EBC
    );
}


void
dgCompressibleMachInletKnudsenOutletCouplingBoundaryField::checkPatchType()
const
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
