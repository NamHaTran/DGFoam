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

#include "rhoBasedConservative.H"
#include "addToRunTimeSelectionTable.H"
#include "eqnOfState.H"
#include "thermoLaw.H"
#include "transportLaw.H"
#include "IOstreams.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rhoBasedConservative, 0);

addToRunTimeSelectionTable(dgThermoConservative, rhoBasedConservative, dictionary);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::rhoBasedConservative::rhoBasedConservative
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgThermoConservative(name, dict, mesh)
{
    // Build sub-models and initialise thermo fields
    initModels();
}


// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

void Foam::rhoBasedConservative::initModels()
{
    // Parse sub-dicts
    const dictionary& dgThermoDict = dict_.subDict("dgThermo");
    const dictionary& mixDict      = dict_.subDict("mixture");

    // --- eqnOfState
    {
        const word eosType(dgThermoDict.lookup("equationOfState"));

        // Create EOS model and transfer ownership into autoPtr
        eqnState_.reset
        (
            eqnOfState::New(eosType, mixDict, mesh_).ptr()
        );
    }

    // --- thermoLaw
    {
        const word thType(dgThermoDict.lookup("thermo"));

        // Create thermo law coupled to EOS model
        thermo_.reset
        (
            thermoLaw::New(thType, mixDict, mesh_, eqnState_).ptr()
        );
    }

    // --- transportLaw (e.g. Sutherland / powerVHS)
    {
        const word trType(dgThermoDict.lookup("transport"));

        // Create transport model coupled to thermo model
        transport_.reset
        (
            transportLaw::New(trType, mixDict, mesh_, thermo_).ptr()
        );
    }

    // --- energy type (e.g. sensible internal energy)
    {
        const word eType(dgThermoDict.lookup("energy"));

        // Create energy model coupled to thermo model
        energy_.reset
        (
            energy::New(eType, mixDict, mesh_, thermo_).ptr()
        );
    }

    // Validate combo
    validateModels();

    // Correct
    // Gas constant from EOS
    {
        R_ = eqnState_().R();

        // Loop over all cells and reconstruct primitive thermo
        for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
        {
            update(cellI);
        }
    }

    // Report selected models
    Info<< "rhoBasedConservative: constructed with models:" << nl
        << "  eqnOfState      = " << eqnState_().type() << nl
        << "  thermo          = " << thermo_().type() << nl
        << "  transport       = " << transport_().type() << nl
        << "  energy          = " << energy_().type() << nl
        << endl;
}


void Foam::rhoBasedConservative::validateModels()
{
    const eqnOfState& eos  = eqnState_();
    const thermoLaw& th    = thermo_();
    const transportLaw& tr = transport_();
    const energy& eng      = energy_();

    // Equation of state: ideal-gas or thermal perfect gas
    if (!(eos.isIdealGas() || eos.isThermalPerfectGas()))
    {
        FatalErrorInFunction
            << "rhoBasedConservative requires an ideal-gas-based or "
            << "thermal-perfect-gas equation of state." << nl
            << "Detected EOS type: " << eos.type()
            << exit(FatalError);
    }

    // Thermo law: perfect-gas or kinetic thermo
    if (!(th.isPerfectGasThermo() || th.isKineticThermo()))
    {
        FatalErrorInFunction
            << "rhoBasedConservative requires a perfect-gas or kinetic "
            << "thermoLaw model." << nl
            << "Detected thermo type: " << th.type()
            << exit(FatalError);
    }

    // Energy: internal-energy-based formulation
    if (!eng.energyInternal())
    {
        FatalErrorInFunction
            << "rhoBasedConservative currently supports only internal-energy "
            << "based energy models (energyInternal() == true)." << nl
            << "Detected energy type: " << eng.type()
            << exit(FatalError);
    }

    // Transport: disallow real-gas transport (not consistent with setup)
    if (tr.isRealGasTransportLaw())
    {
        FatalErrorInFunction
            << "rhoBasedConservative does not support real-gas transport "
            << "models." << nl
            << "Detected transport type: " << tr.type()
            << exit(FatalError);
    }
}

void Foam::rhoBasedConservative::update(const label& cellI)
{
    // Conservative fields
    const GaussField<scalar>& rhoG  = rho_.gaussFields()[cellI];
    const GaussField<vector>& rhoUG = rhoU_.gaussFields()[cellI];
    const GaussField<scalar>& EG    = E_.gaussFields()[cellI];

    // Thermo fields
    GaussField<scalar>& TG      = T_.gaussFields()[cellI];
    GaussField<scalar>& pG      = p_.gaussFields()[cellI];
    GaussField<scalar>& CpG     = Cp_.gaussFields()[cellI];
    GaussField<scalar>& CvG     = Cv_.gaussFields()[cellI];
    GaussField<scalar>& eG      = e_.gaussFields()[cellI];
    GaussField<scalar>& hG      = h_.gaussFields()[cellI];
    GaussField<scalar>& gammaG  = gamma_.gaussFields()[cellI];
    GaussField<scalar>& kappaG  = kappa_.gaussFields()[cellI];
    GaussField<scalar>& muG     = mu_.gaussFields()[cellI];
    GaussField<scalar>& PrG     = Pr_.gaussFields()[cellI];
    GaussField<scalar>& aG      = a_.gaussFields()[cellI];

    // Velocity at Gauss points
    tmp<GaussField<vector>> tUG = rhoUG/rhoG;

    // Kinetic energy
    tmp<GaussField<scalar>> tK = 0.5*magSqr(tUG());

    // Internal energy
    tmp<GaussField<scalar>> te = EG/rhoG - tK();

    eG = te;

    // Temperature
    thermo_().calcT(cellI, eG, TG);

    // Pressure
    eqnState_().calcPFromRhoT(cellI, rhoG, TG, pG);

    // Thermodynamics
    thermo_().calcCp(cellI, TG, CpG);
    thermo_().calcCv(cellI, TG, CvG);
    thermo_().calcGamma(cellI, CpG, CvG, gammaG);
    thermo_().calcH(cellI, TG, hG);
    thermo_().calcSpeedOfSound(cellI, TG, gammaG, aG);

    // Transport
    transport_().calcMu(cellI, TG, muG);
    transport_().calcKappa(cellI, TG, kappaG);
    transport_().calcPr(cellI, TG, PrG);
}


void Foam::rhoBasedConservative::updateBC(const label& cellI)
{
    // Precondition:
    // - the solver has already imposed primitive ghost values on U/T/p,
    // - the matching ghost conservative state (rho, rhoU, E) has already been
    //   rebuilt from those primitive boundary values.
    //
    // This routine must therefore be called after the primary BC update step.
    // It does not re-impose or recompute the primary boundary fields; it only
    // derives the remaining thermo/transport quantities on the boundary state.
    const GaussField<scalar>& rhoG  = rho_.gaussFields()[cellI];
    const GaussField<scalar>& EG    = E_.gaussFields()[cellI];

    GaussField<vector>& UG     = U_.gaussFields()[cellI];
    GaussField<scalar>& TG     = T_.gaussFields()[cellI];
    GaussField<scalar>& pG     = p_.gaussFields()[cellI];
    GaussField<scalar>& CpG    = Cp_.gaussFields()[cellI];
    GaussField<scalar>& CvG    = Cv_.gaussFields()[cellI];
    GaussField<scalar>& eG     = e_.gaussFields()[cellI];
    GaussField<scalar>& hG     = h_.gaussFields()[cellI];
    GaussField<scalar>& gammaG = gamma_.gaussFields()[cellI];
    GaussField<scalar>& kappaG = kappa_.gaussFields()[cellI];
    GaussField<scalar>& muG    = mu_.gaussFields()[cellI];
    GaussField<scalar>& PrG    = Pr_.gaussFields()[cellI];
    GaussField<scalar>& aG     = a_.gaussFields()[cellI];

    // Read back the already-prepared ghost boundary state from the plus side of
    // the face Gauss storage so the derived thermo quantities can be updated.
    const boundaryGaussField<scalar> rhoB =
        rhoG.faceField().extractBCGhostState();
    const boundaryGaussField<scalar> EB =
        EG.faceField().extractBCGhostState();

    boundaryGaussField<vector> UB =
        UG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> TB =
        TG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> pB =
        pG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> CpB =
        CpG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> CvB =
        CvG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> eB =
        eG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> hB =
        hG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> gammaB =
        gammaG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> kappaB =
        kappaG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> muB =
        muG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> PrB =
        PrG.faceField().extractBCGhostState();
    boundaryGaussField<scalar> aB =
        aG.faceField().extractBCGhostState();

    thermo_().calcInternalE(TB, eB);
    thermo_().calcCp(TB, CpB);
    thermo_().calcCv(TB, CvB);
    thermo_().calcGamma(CpB, CvB, gammaB);
    thermo_().calcH(TB, hB);
    thermo_().calcSpeedOfSound(TB, gammaB, aB);

    transport_().calcMu(TB, muB);
    transport_().calcKappa(TB, kappaB);
    transport_().calcPr(TB, PrB);

    // Write back the derived quantities. T and p are assigned unchanged on
    // purpose so the full boundary thermo state stays synchronized in the
    // face storage, but they remain authoritative inputs from the BC stage.
    CpG.faceField().assignBCGhostState(CpB);
    CvG.faceField().assignBCGhostState(CvB);
    eG.faceField().assignBCGhostState(eB);
    hG.faceField().assignBCGhostState(hB);
    gammaG.faceField().assignBCGhostState(gammaB);
    muG.faceField().assignBCGhostState(muB);
    kappaG.faceField().assignBCGhostState(kappaB);
    PrG.faceField().assignBCGhostState(PrB);
    aG.faceField().assignBCGhostState(aB);
}

} // End namespace Foam

// ************************************************************************* //
