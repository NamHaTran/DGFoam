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
#include "eqnOfState.H"   // (expose eqnOfState::New + R())
#include "thermoLaw.H"         // hConst (constant Cp)
#include "transportLaw.H"      // Sutherland, powerVHS
#include "sensibleInternalEnergy.H"
#include "IOstreams.H"
#include <cmath>

namespace Foam
{

// * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rhoBasedConservative, 0);
addToRunTimeSelectionTable(dgThermo, rhoBasedConservative, dictionary);

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

Foam::rhoBasedConservative::rhoBasedConservative
(
    const word& name,
    const dictionary& dict,
    dgGeomMesh& mesh
)
:
    dgThermo(name, dict, mesh),
    Rgas_(Zero)
{
    // 1) Build inner models
    Foam::rhoBasedConservative::initModels();

    // 2) Validate allowed combo (idealGas + hConst + {Sutherland|powerVHS})
    Foam::rhoBasedConservative::validateModelCombo();
}

// * * * * * * * * * * * * * Model Construction  * * * * * * * * * * * * * * //

void Foam::rhoBasedConservative::initModels()
{
    // Expect dictionaries:
    // dgThermo { type rhoBasedConservative; transport sutherland; thermo hConst; eqnOfState idealGas; }
    // mixture  { specie{molWeight ...} thermodynamics{Cp,Hf} eqnOfState{R} transport{As,Ts} }

    // Parse sub-dicts
    const dictionary& dgThermoDict = dict_.subDict("dgThermo");
    const dictionary& mixDict      = dict_.subDict("mixture");

    // --- eqnOfState
    {
        const word eosType(dgThermoDict.lookup("equationOfState"));
        eqnState_.reset( eqnOfState::New(eosType, mixDict).ptr() );
        Rgas_ = eqnState_().R();   // cache
        setR(Rgas_);
    }

    // --- thermoLaw (e.g. hConst = constant Cp)
    {
        const word thType(dgThermoDict.lookup("thermo"));
        thermo_.reset( thermoLaw::New(thType, mixDict).ptr() );
        setCp( thermo_().Cp(300.0) );  // if hConst, Cp is constant; 300 K placeholder
    }

    // --- transportLaw (e.g. Sutherland / powerVHS)
    {
        const word trType(dgThermoDict.lookup("transport"));
        transport_.reset( transportLaw::New(trType, mixDict).ptr() );
    }

    // --- energy type (e.g. sensible internal energy)
    {
        const word eType(dgThermoDict.lookup("energy"));
        energy_.reset( energy::New(eType, mixDict).ptr() );
    }
}

void Foam::rhoBasedConservative::validateModelCombo()
{
    // Allow only the listed inner models (extend as needed)
    const word eosName = eqnState_().type();
    const word thName  = thermo_().type();
    const word trName  = transport_().type();
    const word eName   = energy_().type();

    const bool eosOK = (eosName == "idealGas");
    const bool thOK  = (thName  == "hConst" || thName == "constantCp"); // alias
    const bool trOK  = (trName  == "Sutherland" || trName == "powerVHS");
    const bool eOK   = (eName   == "sensibleInternalEnergy");

    if (!(eosOK && thOK && trOK && eOK))
    {
        FatalErrorInFunction
            << "rhoBasedConservative: unsupported inner model combo:" << nl
            << "  eqnOfState      = " << eosName << nl
            << "  thermo          = " << thName  << nl
            << "  transport       = " << trName  << nl
            << "  energy          = " << eName   << nl
            << "Allowed: idealGas + hConst(constantCp) + {Sutherland|powerVHS} + sensibleInternalEnergy"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * Pipeline (update)  * * * * * * * * * * * * * * * //

void Foam::rhoBasedConservative::update
(
    const dgThermoInputs& in,
    dgThermoOutputs& out
)
{
    // 0) Sanity: required inputs
    if (!in.rhoC || !in.rhoU || !in.rhoE)
    {
        FatalErrorInFunction
            << "rhoBasedConservative requires conserved inputs: rhoC, rhoU, rhoE." << nl
            << exit(FatalError);
    }

    // 1) Recover primitive from conserved: rho, U, e
    scalar rho = max(*in.rhoC, SMALL);
    vector U   = (*in.rhoU)/rho;

    const scalar E   = (*in.rhoE)/rho;
    const scalar kin = 0.5*magSqr(U);
    scalar e         = E - kin;

    // 2) Get thermo constant
    const scalar Cp = thermo_().Cp(300.0);      // hConst -> const

    // 3) Get T
    scalar T = energy_().T(e, Cp - Rgas_, Cp );

    // 4) Thermo: Cp, gamma  (hConst -> Cp = const; if not const, add T-iteration)
    
    const scalar gamma = Cp/(Cp - Rgas_);
    const scalar h = thermo_().h(T);

    // 5) EOS (ideal gas energy form): p, T, a
    const scalar p = eqnState_().p(rho,T);
    const scalar a = std::sqrt(gamma*p/rho);

    // 6) Transport: mu(T) and kappa(Cp,T) via base transportLaw
    const scalar mu    = transport_().mu(T);
    const scalar kappa = transport_().kappa(Cp, T);
    const scalar Pr    = transport_().Pr(T);


    // 5) Cache to base (optional)
    setCp(Cp); setE(e); setH(h); setMu(mu); setKappa(kappa); setPr(Pr); setA(a); setGamma(gamma);

    // 6) Write outputs if pointers are provided
    if (out.rho)   *out.rho   = rho;
    if (out.U)     *out.U     = U;
    if (out.p)     *out.p     = p;
    if (out.T)     *out.T     = T;

    if (out.a)     *out.a     = a;
    if (out.mu)    *out.mu    = mu;
    if (out.kappa) *out.kappa = kappa;
    if (out.Cp)    *out.Cp    = Cp;
    if (out.h)     *out.h     = h;
    if (out.e)     *out.e     = e;
    if (out.Pr)    *out.Pr    = Pr;
    if (out.gamma) *out.gamma = gamma;
}

// ************************************************************************* //

} // End namespace Foam
