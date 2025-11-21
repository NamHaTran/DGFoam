/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
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
    // Conservative fields at Gauss points
    const GaussField<scalar>& rhoG  = rho_.gaussFields()[cellI];
    const GaussField<vector>& rhoUG = rhoU_.gaussFields()[cellI];
    const GaussField<scalar>& EG    = E_.gaussFields()[cellI];

    // Thermo fields at Gauss points (references into storage)
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

    // Reconstruct velocity at Gauss points: U = rhoU / rho
    GaussField<vector> UG = rhoUG/rhoG;

    // Kinetic energy per unit mass at Gauss points
    GaussField<scalar> kineticE = 0.5*magSqr(UG);

    // Internal energy per unit mass:
    //   e = E/rho - 0.5|U|^2
    GaussField<scalar> internalE = EG/rhoG - kineticE;
    eG = internalE;

    // Temperature from internal energy
    thermo_().calcT(cellI, internalE, TG);

    // Pressure from EOS: p = p(rho, T)
    eqnState_().calcPFromRhoT(cellI, rhoG, TG, pG);

    // Thermodynamic properties
    thermo_().calcCp(cellI, TG, CpG);
    thermo_().calcCv(cellI, TG, CvG);
    thermo_().calcGamma(cellI, CpG, CvG, gammaG);
    thermo_().calcH(cellI, TG, hG);
    thermo_().calcSpeedOfSound(cellI, TG, gammaG, aG);

    // Transport properties
    transport_().calcMu(cellI, TG, muG);
    transport_().calcKappa(cellI, TG, kappaG);
    transport_().calcPr(cellI, TG, PrG);
}

} // End namespace Foam

// ************************************************************************* //