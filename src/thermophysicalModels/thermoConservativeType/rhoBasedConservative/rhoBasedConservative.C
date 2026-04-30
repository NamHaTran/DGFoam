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
#include "dgCompressibleBoundaryManager.H"

namespace Foam
{
namespace
{
struct RhoHeState
{
    scalar rho;
    scalar he;
};


struct RhoTState
{
    scalar rho;
    scalar T;
};


struct GradRhoHeState
{
    vector gradRho;
    vector gradHe;
};


struct GradRhoTState
{
    vector gradRho;
    vector gradT;
};


scalar safeDensity(const scalar rho)
{
    return mag(rho) > VSMALL ? rho : (rho >= 0 ? VSMALL : -VSMALL);
}


vector gradKineticEnergy
(
    const vector& U,
    const tensor& gradU
)
{
    return vector
    (
        U.x()*gradU.xx() + U.y()*gradU.yx() + U.z()*gradU.zx(),
        U.x()*gradU.xy() + U.y()*gradU.yy() + U.z()*gradU.zy(),
        U.x()*gradU.xz() + U.y()*gradU.yz() + U.z()*gradU.zz()
    );
}
}

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
        << "  selfDiffusion   = " << Switch(selfDiffusion_) << nl
        << "  selfDiffusionD0 = " << selfDiffCoeffScale_ << nl
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
    const GaussField<scalar>& EG    = E_.gaussFields()[cellI];
    const GaussField<vector>& UG    = U_.gaussFields()[cellI];

    // Thermo fields
    GaussField<scalar>& TG      = T_.gaussFields()[cellI];
    GaussField<scalar>& pG      = p_.gaussFields()[cellI];
    GaussField<scalar>& heG     = he_.gaussFields()[cellI];
    GaussField<scalar>& aG      = a_.gaussFields()[cellI];
    GaussField<scalar>& muG     = mu_.gaussFields()[cellI];
    GaussField<scalar>& kappaG  = kappa_.gaussFields()[cellI];

    const auto rhoExpr = dg::expr(rhoG);
    const auto heExpr =
        dg::expr(EG)/rhoExpr - 0.5*dg::magSqr(dg::expr(UG));

    dg::assign
    (
        heG,
        heExpr
    );

    // Rebuild thermo state from the current conservative polynomial. The
    // expression assignments write cell and face Gauss storage in one path and
    // avoid materialising intermediate fields.
    dg::assign
    (
        TG,
        dg::map
        (
            [this](const scalar rho, const scalar he)
            {
                return clampTemperature
                (
                    calcRawTemperatureFromRhoHe(rho, he),
                    scalar(SMALL)
                );
            },
            rhoExpr,
            dg::expr(heG)
        )
    );

    dg::assign
    (
        pG,
        dg::map
        (
            [this](const scalar rho, const scalar T)
            {
                return eos().calcPFromRhoT(rho, T);
            },
            rhoExpr,
            dg::expr(TG)
        )
    );

    dg::assign
    (
        aG,
        dg::map
        (
            [this](const scalar T)
            {
                const scalar Cp = thermo().calcCp(T);
                const scalar Cv = thermo().calcCv(T);
                return thermo().calcSpeedOfSound
                (
                    T,
                    thermo().calcGamma(Cp, Cv)
                );
            },
            dg::expr(TG)
        )
    );

    dg::assign
    (
        muG,
        dg::map
        (
            [this](const scalar T)
            {
                return transport().calcMu(T);
            },
            dg::expr(TG)
        )
    );

    dg::assign
    (
        kappaG,
        dg::map
        (
            [this](const scalar mu, const scalar T)
            {
                return mu*thermo().calcCp(T)/transport().calcPr(T);
            },
            dg::expr(muG),
            dg::expr(TG)
        )
    );
}


void Foam::rhoBasedConservative::updateGradient(const label& cellI)
{
    const GaussField<scalar>& rhoG = rho_.gaussFields()[cellI];
    const GaussField<scalar>& EG = E_.gaussFields()[cellI];
    const GaussField<vector>& UG = U_.gaussFields()[cellI];
    const GaussField<scalar>& TG = T_.gaussFields()[cellI];
    const GaussField<scalar>& heG = he_.gaussFields()[cellI];
    const GaussField<vector>& SRhoG = SRho_.gaussFields()[cellI];
    const GaussField<vector>& SEG = SE_.gaussFields()[cellI];
    const GaussField<tensor>& gradUG = gradU_.gaussFields()[cellI];

    GaussField<vector>& gradPG = gradP_.gaussFields()[cellI];
    GaussField<vector>& gradTG = gradT_.gaussFields()[cellI];

    const auto rhoExpr = dg::expr(rhoG);
    const auto rhoSafeExpr =
        dg::map
        (
            [](const scalar rho)
            {
                return safeDensity(rho);
            },
            rhoExpr
        );
    const auto gradRhoExpr = dg::expr(SRhoG);
    const auto gradHeExpr =
        dg::expr(SEG)/rhoSafeExpr
      - (dg::expr(EG)/(rhoSafeExpr*rhoSafeExpr))*gradRhoExpr
      - dg::map
        (
            [](const vector& U, const tensor& gradU)
            {
                return gradKineticEnergy(U, gradU);
            },
            dg::expr(UG),
            dg::expr(gradUG)
        );

    const auto rhoHeExpr =
        dg::map
        (
            [](const scalar rho, const scalar he)
            {
                return RhoHeState{rho, he};
            },
            rhoExpr,
            dg::expr(heG)
        );

    const auto gradRhoHeExpr =
        dg::map
        (
            [](const vector& gradRho, const vector& gradHe)
            {
                return GradRhoHeState{gradRho, gradHe};
            },
            gradRhoExpr,
            gradHeExpr
        );

    const bool useEosGradT = heIsInternalEnergy() && eos().canCalcTFromRhoE();

    dg::assign
    (
        gradTG,
        dg::map
        (
            [this, useEosGradT]
            (
                const RhoHeState& state,
                const GradRhoHeState& gradState
            )
            {
                if (useEosGradT)
                {
                    return eos().calcGradTFromRhoE
                    (
                        state.rho,
                        state.he,
                        gradState.gradRho,
                        gradState.gradHe
                    );
                }

                return energyModel().calcGradTfromHe
                (
                    state.he,
                    gradState.gradHe
                );
            },
            rhoHeExpr,
            gradRhoHeExpr
        )
    );

    const auto rhoTExpr =
        dg::map
        (
            [](const scalar rho, const scalar T)
            {
                return RhoTState{rho, T};
            },
            rhoExpr,
            dg::expr(TG)
        );

    const auto gradRhoTExpr =
        dg::map
        (
            [](const vector& gradRho, const vector& gradT)
            {
                return GradRhoTState{gradRho, gradT};
            },
            gradRhoExpr,
            dg::expr(gradTG)
        );

    dg::assign
    (
        gradPG,
        dg::map
        (
            [this](const RhoTState& state, const GradRhoTState& gradState)
            {
                return eos().calcGradPFromRhoT
                (
                    state.rho,
                    state.T,
                    gradState.gradRho,
                    gradState.gradT
                );
            },
            rhoTExpr,
            gradRhoTExpr
        )
    );

    // Correct thermo-derived gradients on boundary faces using the boundary
    // manager so wall/outlet BCs can enforce their face-normal constraints on
    // both grad(T) and grad(p).
    const objectRegistry& obr = mesh_.getFvMesh();

    if (obr.foundObject<dgCompressibleBoundaryManager>("dgCompressibleBoundaryManager"))
    {
        const dgCompressibleBoundaryManager& compressibleBC =
            obr.lookupObject<dgCompressibleBoundaryManager>
            (
                "dgCompressibleBoundaryManager"
            );

        compressibleBC.correctGradient(gradTG);
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
    GaussField<scalar>& muG     = mu_.gaussFields()[cellI];
    GaussField<scalar>& kappaG  = kappa_.gaussFields()[cellI];

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
    boundaryGaussField<scalar> muB(rhoB.size());
    boundaryGaussField<scalar> kappaB(rhoB.size());

    // Rebuild thermo boundary traces from the ghost conservative state.
    for (label gpI = 0; gpI < rhoB.size(); ++gpI)
    {
        const vector U = rhoUB[gpI]/rhoB[gpI];
        const scalar k = 0.5*magSqr(U);

        heB[gpI] = EB[gpI]/rhoB[gpI] - k;

        const scalar T =
            clampTemperature
            (
                calcRawTemperatureFromRhoHe(rhoB[gpI], heB[gpI]),
                scalar(SMALL)
            );

        TB[gpI] = T;
        pB[gpI] = eos().calcPFromRhoT(rhoB[gpI], T);

        const scalar Cp = thermo().calcCp(T);
        const scalar Cv = thermo().calcCv(T);
        const scalar gamma = thermo().calcGamma(Cp, Cv);
        aB[gpI] = thermo().calcSpeedOfSound(T, gamma);
        muB[gpI] = transport().calcMu(T);
        kappaB[gpI] = muB[gpI]*Cp/transport().calcPr(T);
    }

    TG.faceField().assignBCGhostState(TB);
    pG.faceField().assignBCGhostState(pB);
    heG.faceField().assignBCGhostState(heB);
    aG.faceField().assignBCGhostState(aB);
    muG.faceField().assignBCGhostState(muB);
    kappaG.faceField().assignBCGhostState(kappaB);
}

} // End namespace Foam

// ************************************************************************* //
