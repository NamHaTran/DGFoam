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
#include "dgExpr.H"

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

        // Populate plus-side processor traces for thermo-owned scalar fields
        // so the very first stage starts from a consistent parallel state.
        synch();
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

    // Equation of state: ideal-gas, thermal-perfect gas, or supported real gas
    if (!(eos.isIdealGas() || eos.isThermalPerfectGas() || eos.isRealGas()))
    {
        FatalErrorInFunction
            << "rhoBasedConservative requires an ideal-gas, thermal-perfect, "
            << "or supported real-gas equation of state." << nl
            << "Detected EOS type: " << eos.type()
            << exit(FatalError);
    }

    // Real-gas support currently uses constantCp only as a reference caloric
    // relation; ideal-gas and kinetic paths retain the previous flexibility.
    if (eos.isRealGas())
    {
        if (th.type() != "constantCp")
        {
            FatalErrorInFunction
                << "rhoBasedConservative currently supports real-gas EOS only "
                << "with thermoLaw 'constantCp'." << nl
                << "Detected thermo type: " << th.type()
                << exit(FatalError);
        }

        if
        (
            !eos.canCalcTFromRhoE()
         || !eos.canCalcPFromRhoE()
         || !eos.canCalcAFromRhoE()
         || !eos.canCalcEFromRhoT()
        )
        {
            FatalErrorInFunction
                << "Selected real-gas EOS '" << eos.type()
                << "' does not provide the conservative reconstruction hooks "
                << "required by rhoBasedConservative."
                << exit(FatalError);
        }
    }
    else if (!(th.isPerfectGasThermo() || th.isKineticThermo()))
    {
        FatalErrorInFunction
            << "rhoBasedConservative requires a perfect-gas or kinetic "
            << "thermoLaw model." << nl
            << "Detected thermo type: " << th.type()
            << exit(FatalError);
    }

    // Energy: currently only internal-energy-based formulation is supported
    if (!eng.heIsInternalEnergy())
    {
        FatalErrorInFunction
            << "rhoBasedConservative currently supports only internal-energy "
            << "based energy models (heType == internalEnergy)." << nl
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
    GaussField<scalar>& heG     = he_.gaussFields()[cellI];
    GaussField<scalar>& aG      = a_.gaussFields()[cellI];

    const auto velocityExpr = dg::expr(rhoUG)/dg::expr(rhoG);

    dg::assign
    (
        heG,
        dg::expr(EG)/dg::expr(rhoG) - 0.5*dg::magSqr(velocityExpr)
    );

    cellGaussField<scalar>& aCell = aG.cellField();
    cellGaussField<scalar>& TCellMutable = TG.cellField();
    cellGaussField<scalar>& pCellMutable = pG.cellField();
    const cellGaussField<scalar>& rhoCell = rhoG.cellField();
    const cellGaussField<scalar>& heCell = heG.cellField();

    for (label gpI = 0; gpI < heCell.size(); ++gpI)
    {
        const scalar rho = rhoCell[gpI];
        const scalar he = heCell[gpI];
        const scalar T = calcTemperatureFromRhoHe(rho, he);

        TCellMutable[gpI] = T;
        pCellMutable[gpI] = calcPressureFromRhoHe(rho, he);
        aCell[gpI] = calcSpeedOfSoundFromRhoHe(rho, he);
    }

    faceGaussField<scalar>& aFace = aG.faceField();
    faceGaussField<scalar>& TFaceMutable = TG.faceField();
    faceGaussField<scalar>& pFaceMutable = pG.faceField();
    const faceGaussField<scalar>& rhoFace = rhoG.faceField();
    const faceGaussField<scalar>& heFace = heG.faceField();

    for (label gpI = 0; gpI < heFace.nGauss(); ++gpI)
    {
        const scalar rhoMinus = rhoFace.minusValue(gpI);
        const scalar heMinus = heFace.minusValue(gpI);
        const scalar TMinus = calcTemperatureFromRhoHe(rhoMinus, heMinus);

        TFaceMutable.minusValueAt(gpI) = TMinus;
        pFaceMutable.minusValueAt(gpI) =
            calcPressureFromRhoHe(rhoMinus, heMinus);
        aFace.minusValueAt(gpI) =
            calcSpeedOfSoundFromRhoHe(rhoMinus, heMinus);

        const scalar rhoPlus = rhoFace.plusValue(gpI);
        const scalar hePlus = heFace.plusValue(gpI);
        const scalar TPlus = calcTemperatureFromRhoHe(rhoPlus, hePlus);

        TFaceMutable.plusValueAt(gpI) = TPlus;
        pFaceMutable.plusValueAt(gpI) =
            calcPressureFromRhoHe(rhoPlus, hePlus);
        aFace.plusValueAt(gpI) =
            calcSpeedOfSoundFromRhoHe(rhoPlus, hePlus);
    }
}


void Foam::rhoBasedConservative::updateBC(const label& cellI)
{
    // The boundary manager has already written the authoritative ghost state
    // into rho/rhoU/E. Rebuild the thermal boundary traces from that
    // conservative data so flux assembly sees a fully consistent plus state.
    const GaussField<scalar>& rhoG  = rho_.gaussFields()[cellI];
    const GaussField<vector>& rhoUG = rhoU_.gaussFields()[cellI];
    const GaussField<scalar>& EG    = E_.gaussFields()[cellI];

    GaussField<scalar>& TG      = T_.gaussFields()[cellI];
    GaussField<scalar>& pG      = p_.gaussFields()[cellI];
    GaussField<scalar>& heG     = he_.gaussFields()[cellI];
    GaussField<scalar>& aG      = a_.gaussFields()[cellI];

    const boundaryGaussField<scalar> rhoB =
        rhoG.faceField().extractBCGhostState();
    const boundaryGaussField<vector> rhoUB =
        rhoUG.faceField().extractBCGhostState();
    const boundaryGaussField<scalar> EB =
        EG.faceField().extractBCGhostState();

    boundaryGaussField<scalar> TB(rhoB.size());
    boundaryGaussField<scalar> pB(rhoB.size());
    boundaryGaussField<scalar> heB(rhoB.size());
    boundaryGaussField<scalar> aB(rhoB.size());

    for (label gpI = 0; gpI < rhoB.size(); ++gpI)
    {
        const vector U = rhoUB[gpI]/rhoB[gpI];
        const scalar k = 0.5*magSqr(U);

        heB[gpI] = EB[gpI]/rhoB[gpI] - k;
        TB[gpI] = calcTemperatureFromRhoHe(rhoB[gpI], heB[gpI]);
        pB[gpI] = calcPressureFromRhoHe(rhoB[gpI], heB[gpI]);
        aB[gpI] = calcSpeedOfSoundFromRhoHe(rhoB[gpI], heB[gpI]);
    }

    TG.faceField().assignBCGhostState(TB);
    pG.faceField().assignBCGhostState(pB);
    heG.faceField().assignBCGhostState(heB);
    aG.faceField().assignBCGhostState(aB);
}

} // End namespace Foam

// ************************************************************************* //
