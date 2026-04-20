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
namespace
{
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


vector gradSpecificInternalEnergy
(
    const scalar rho,
    const vector& U,
    const scalar E,
    const vector& gradRho,
    const tensor& gradU,
    const vector& gradE
)
{
    const scalar rhoSafe = safeDensity(rho);

    return
        gradE/rhoSafe
      - (E/sqr(rhoSafe))*gradRho
      - gradKineticEnergy(U, gradU);
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


scalar Foam::rhoBasedConservative::calcTMean
(
    const scalar rho0,
    const vector& rhoU0,
    const scalar E0
) const
{
    const scalar rho0Safe =
        mag(rho0) > VSMALL ? rho0 : (rho0 >= 0 ? VSMALL : -VSMALL);
    const vector U0 = rhoU0/rho0Safe;
    const scalar he0 = E0/rho0Safe - 0.5*magSqr(U0);

    return clampTemperature(calcRawTemperatureFromRhoHe(rho0Safe, he0), scalar(SMALL));
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
    GaussField<scalar>& muG     = mu_.gaussFields()[cellI];

    const auto velocityExpr = dg::expr(rhoUG)/dg::expr(rhoG);

    dg::assign
    (
        heG,
        dg::expr(EG)/dg::expr(rhoG) - 0.5*dg::magSqr(velocityExpr)
    );

    // Reconstruct one cell-mean temperature from the conservative mean
    // state (DoF 0). This acts as a fallback when a pointwise reconstruction
    // produces a slightly negative temperature during restart/rebuild.
    const scalar TMean =
        calcTMean
        (
            rho_.dof()[cellI][0],
            rhoU_.dof()[cellI][0],
            E_.dof()[cellI][0]
        );
    setTMean(cellI, TMean);

    cellGaussField<scalar>& aCell = aG.cellField();
    cellGaussField<scalar>& muCell = muG.cellField();
    cellGaussField<scalar>& TCellMutable = TG.cellField();
    cellGaussField<scalar>& pCellMutable = pG.cellField();
    const cellGaussField<scalar>& rhoCell = rhoG.cellField();
    const cellGaussField<scalar>& heCell = heG.cellField();

    // Rebuild T, p, and a at cell-interior Gauss points from the current
    // conservative polynomial. If T becomes negative, fall back to TMean.
    for (label gpI = 0; gpI < heCell.size(); ++gpI)
    {
        const scalar rho = rhoCell[gpI];
        const scalar he = heCell[gpI];
        const scalar T = calcTemperatureFromRhoHe(cellI, rho, he);

        TCellMutable[gpI] = T;
        pCellMutable[gpI] = eos().calcPFromRhoT(rho, T);

        const scalar Cp = thermo().calcCp(T);
        const scalar Cv = thermo().calcCv(T);
        const scalar gamma = thermo().calcGamma(Cp, Cv);
        aCell[gpI] = thermo().calcSpeedOfSound(T, gamma);
        muCell[gpI] = transport().calcMu(T);
    }

    faceGaussField<scalar>& aFace = aG.faceField();
    faceGaussField<scalar>& muFace = muG.faceField();
    faceGaussField<scalar>& TFaceMutable = TG.faceField();
    faceGaussField<scalar>& pFaceMutable = pG.faceField();
    const faceGaussField<scalar>& rhoFace = rhoG.faceField();
    const faceGaussField<scalar>& heFace = heG.faceField();

    // Rebuild face thermo states on both minus and plus sides. For the plus
    // side of internal faces, also compute the neighbour-cell mean
    // temperature from the neighbour DoF-0 state as its local fallback.
    for (label localFaceI = 0; localFaceI < heFace.nFaces(); ++localFaceI)
    {
        scalar TMeanPlus = TMean;

        if (!heFace.isBoundary(localFaceI))
        {
            TMeanPlus =
                calcTMean
                (
                    rho_.neighbourCellDof(cellI, localFaceI)[0],
                    rhoU_.neighbourCellDof(cellI, localFaceI)[0],
                    E_.neighbourCellDof(cellI, localFaceI)[0]
                );
        }

        // Traverse the Gauss points on the current face and rebuild minus/plus
        // thermo traces from conservative traces. Small negative temperatures
        // are replaced by the corresponding mean-state temperature.
        for (label gpI = 0; gpI < heFace.nGaussPerFace(); ++gpI)
        {
            const scalar rhoMinus = rhoFace.minusValueOnFace(localFaceI, gpI);
            const scalar heMinus = heFace.minusValueOnFace(localFaceI, gpI);
            const scalar TMinus =
                clampTemperature
                (
                    calcRawTemperatureFromRhoHe(rhoMinus, heMinus),
                    TMean
                );

            TFaceMutable.minusValueOnFace(localFaceI, gpI) = TMinus;
            pFaceMutable.minusValueOnFace(localFaceI, gpI) =
                eos().calcPFromRhoT(rhoMinus, TMinus);

            {
                const scalar Cp = thermo().calcCp(TMinus);
                const scalar Cv = thermo().calcCv(TMinus);
                const scalar gamma = thermo().calcGamma(Cp, Cv);
                aFace.minusValueOnFace(localFaceI, gpI) =
                    thermo().calcSpeedOfSound(TMinus, gamma);
                muFace.minusValueOnFace(localFaceI, gpI) =
                    transport().calcMu(TMinus);
            }

            const scalar rhoPlus = rhoFace.plusValueOnFace(localFaceI, gpI);
            const scalar hePlus = heFace.plusValueOnFace(localFaceI, gpI);
            const scalar TPlus =
                clampTemperature
                (
                    calcRawTemperatureFromRhoHe(rhoPlus, hePlus),
                    TMeanPlus
                );

            TFaceMutable.plusValueOnFace(localFaceI, gpI) = TPlus;
            pFaceMutable.plusValueOnFace(localFaceI, gpI) =
                eos().calcPFromRhoT(rhoPlus, TPlus);

            {
                const scalar Cp = thermo().calcCp(TPlus);
                const scalar Cv = thermo().calcCv(TPlus);
                const scalar gamma = thermo().calcGamma(Cp, Cv);
                aFace.plusValueOnFace(localFaceI, gpI) =
                    thermo().calcSpeedOfSound(TPlus, gamma);
                muFace.plusValueOnFace(localFaceI, gpI) =
                    transport().calcMu(TPlus);
            }
        }
    }
}


void Foam::rhoBasedConservative::updateGradient(const label& cellI)
{
    // Step 1: gather the already-reconstructed primitive/thermal state.
    // rho/E are conservative inputs, while U/T/he were refreshed earlier by
    // update(). This function assumes that sequence has already happened.
    const GaussField<scalar>& rhoG = rho_.gaussFields()[cellI];
    const GaussField<scalar>& EG = E_.gaussFields()[cellI];
    const GaussField<vector>& UG = U_.gaussFields()[cellI];
    const GaussField<scalar>& TG = T_.gaussFields()[cellI];
    const GaussField<scalar>& heG = he_.gaussFields()[cellI];

    // Step 2: gather the gradient inputs. The conservative gradients
    // SRho/SE come from the conservative solution, and gradU is intentionally
    // computed outside thermo so the thermo layer only applies thermodynamic
    // chain rules.
    const GaussField<vector>& SRhoG = SRho_.gaussFields()[cellI];
    const GaussField<vector>& SEG = SE_.gaussFields()[cellI];
    const GaussField<tensor>& gradUG = gradU_.gaussFields()[cellI];

    // Step 3: get writable storage for the thermo-derived gradients.
    // updateGradient() owns only gradP and gradT; it does not modify gradU or
    // the conservative gradients.
    GaussField<vector>& gradPG = gradP_.gaussFields()[cellI];
    GaussField<vector>& gradTG = gradT_.gaussFields()[cellI];

    // Step 4: work on interior Gauss points first. These views avoid repeated
    // cellField() lookups in the point loop below.
    const cellGaussField<scalar>& rhoCell = rhoG.cellField();
    const cellGaussField<scalar>& ECell = EG.cellField();
    const cellGaussField<vector>& UCell = UG.cellField();
    const cellGaussField<scalar>& TCell = TG.cellField();
    const cellGaussField<scalar>& heCell = heG.cellField();
    const cellGaussField<vector>& SRhoCell = SRhoG.cellField();
    const cellGaussField<vector>& SECell = SEG.cellField();
    const cellGaussField<tensor>& gradUCell = gradUG.cellField();

    cellGaussField<vector>& gradPCell = gradPG.cellField();
    cellGaussField<vector>& gradTCell = gradTG.cellField();

    // Step 5: reconstruct gradients at cell-interior Gauss points.
    //
    // For rhoBasedConservative the stored thermal energy is the specific
    // internal energy
    //
    //     he = E/rho - 0.5*magSqr(U)
    //
    // so its gradient is obtained from conservative gradients plus the
    // externally supplied velocity gradient. Once grad(he) is known, the
    // selected EOS/energy model maps it to grad(T), then EOS maps
    // (rho, T, grad(rho), grad(T)) to grad(p).
    for (label gpI = 0; gpI < heCell.size(); ++gpI)
    {
        // Step 5a: apply the conservative-to-internal-energy chain rule.
        const vector gradHe =
            gradSpecificInternalEnergy
            (
                rhoCell[gpI],
                UCell[gpI],
                ECell[gpI],
                SRhoCell[gpI],
                gradUCell[gpI],
                SECell[gpI]
            );

        if (heIsInternalEnergy() && eos().canCalcTFromRhoE())
        {
            // Step 5b-1: real-gas or EOS-owned caloric closures can provide
            // grad(T) directly from T(rho, he).
            gradTCell[gpI] =
                eos().calcGradTFromRhoE
                (
                    rhoCell[gpI],
                    heCell[gpI],
                    SRhoCell[gpI],
                    gradHe
                );
        }
        else
        {
            // Step 5b-2: otherwise use the energy model relation T(he).
            // The base energy API falls back to finite differences when a
            // concrete model has no analytic implementation.
            gradTCell[gpI] =
                energyModel().calcGradTfromHe(heCell[gpI], gradHe);
        }

        // Step 5c: pressure is always differentiated through p(rho, T).
        // EOS models with explicit derivatives override this hook; the base
        // EOS implementation provides a finite-difference fallback.
        gradPCell[gpI] =
            eos().calcGradPFromRhoT
            (
                rhoCell[gpI],
                TCell[gpI],
                SRhoCell[gpI],
                gradTCell[gpI]
            );
    }

    // Step 6: repeat the same chain-rule reconstruction on face Gauss data.
    // Both minus and plus traces are updated because viscous/gradient-aware
    // fluxes may read either side of an internal or processor face.
    const faceGaussField<scalar>& rhoFace = rhoG.faceField();
    const faceGaussField<scalar>& EFace = EG.faceField();
    const faceGaussField<vector>& UFace = UG.faceField();
    const faceGaussField<scalar>& TFace = TG.faceField();
    const faceGaussField<scalar>& heFace = heG.faceField();
    const faceGaussField<vector>& SRhoFace = SRhoG.faceField();
    const faceGaussField<vector>& SEFace = SEG.faceField();
    const faceGaussField<tensor>& gradUFace = gradUG.faceField();

    faceGaussField<vector>& gradPFace = gradPG.faceField();
    faceGaussField<vector>& gradTFace = gradTG.faceField();

    for (label localFaceI = 0; localFaceI < heFace.nFaces(); ++localFaceI)
    {
        for (label gpI = 0; gpI < heFace.nGaussPerFace(); ++gpI)
        {
            // Step 6a: minus-side trace owned by the current cell.
            const vector gradHeMinus =
                gradSpecificInternalEnergy
                (
                    rhoFace.minusValueOnFace(localFaceI, gpI),
                    UFace.minusValueOnFace(localFaceI, gpI),
                    EFace.minusValueOnFace(localFaceI, gpI),
                    SRhoFace.minusValueOnFace(localFaceI, gpI),
                    gradUFace.minusValueOnFace(localFaceI, gpI),
                    SEFace.minusValueOnFace(localFaceI, gpI)
                );

            if (heIsInternalEnergy() && eos().canCalcTFromRhoE())
            {
                // Step 6b: convert minus-side grad(he) to grad(T).
                gradTFace.minusValueOnFace(localFaceI, gpI) =
                    eos().calcGradTFromRhoE
                    (
                        rhoFace.minusValueOnFace(localFaceI, gpI),
                        heFace.minusValueOnFace(localFaceI, gpI),
                        SRhoFace.minusValueOnFace(localFaceI, gpI),
                        gradHeMinus
                    );
            }
            else
            {
                // Step 6b fallback path for models where T is supplied by
                // the energy model instead of the EOS.
                gradTFace.minusValueOnFace(localFaceI, gpI) =
                    energyModel().calcGradTfromHe
                    (
                        heFace.minusValueOnFace(localFaceI, gpI),
                        gradHeMinus
                    );
            }

            // Step 6c: convert minus-side grad(rho), grad(T) to grad(p).
            gradPFace.minusValueOnFace(localFaceI, gpI) =
                eos().calcGradPFromRhoT
                (
                    rhoFace.minusValueOnFace(localFaceI, gpI),
                    TFace.minusValueOnFace(localFaceI, gpI),
                    SRhoFace.minusValueOnFace(localFaceI, gpI),
                    gradTFace.minusValueOnFace(localFaceI, gpI)
                );

            // Step 6d: plus-side trace. On internal/processor faces this is
            // the neighbour trace already stored in the faceGaussField.
            const vector gradHePlus =
                gradSpecificInternalEnergy
                (
                    rhoFace.plusValueOnFace(localFaceI, gpI),
                    UFace.plusValueOnFace(localFaceI, gpI),
                    EFace.plusValueOnFace(localFaceI, gpI),
                    SRhoFace.plusValueOnFace(localFaceI, gpI),
                    gradUFace.plusValueOnFace(localFaceI, gpI),
                    SEFace.plusValueOnFace(localFaceI, gpI)
                );

            if (heIsInternalEnergy() && eos().canCalcTFromRhoE())
            {
                // Step 6e: convert plus-side grad(he) to grad(T).
                gradTFace.plusValueOnFace(localFaceI, gpI) =
                    eos().calcGradTFromRhoE
                    (
                        rhoFace.plusValueOnFace(localFaceI, gpI),
                        heFace.plusValueOnFace(localFaceI, gpI),
                        SRhoFace.plusValueOnFace(localFaceI, gpI),
                        gradHePlus
                    );
            }
            else
            {
                // Step 6e fallback path through the energy model.
                gradTFace.plusValueOnFace(localFaceI, gpI) =
                    energyModel().calcGradTfromHe
                    (
                        heFace.plusValueOnFace(localFaceI, gpI),
                        gradHePlus
                    );
            }

            // Step 6f: convert plus-side grad(rho), grad(T) to grad(p).
            gradPFace.plusValueOnFace(localFaceI, gpI) =
                eos().calcGradPFromRhoT
                (
                    rhoFace.plusValueOnFace(localFaceI, gpI),
                    TFace.plusValueOnFace(localFaceI, gpI),
                    SRhoFace.plusValueOnFace(localFaceI, gpI),
                    gradTFace.plusValueOnFace(localFaceI, gpI)
                );
        }
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

    const boundaryGaussField<scalar> rhoB =
        rhoG.faceField().extractBCGhostState();
    const boundaryGaussField<vector> rhoUB =
        rhoUG.faceField().extractBCGhostState();
    const boundaryGaussField<scalar> EB =
        EG.faceField().extractBCGhostState();

    // Reconstruct one owner-cell mean temperature from DoF 0 and use it as
    // the boundary fallback when a ghost-state reconstruction yields a small
    // negative temperature.
    const scalar TMean = TMean_[cellI];

    boundaryGaussField<scalar> TB(rhoB.size());
    boundaryGaussField<scalar> pB(rhoB.size());
    boundaryGaussField<scalar> heB(rhoB.size());
    boundaryGaussField<scalar> aB(rhoB.size());
    boundaryGaussField<scalar> muB(rhoB.size());

    // Rebuild thermo boundary traces from the ghost conservative state and
    // replace any negative temperature with the owner-cell mean temperature.
    for (label gpI = 0; gpI < rhoB.size(); ++gpI)
    {
        const vector U = rhoUB[gpI]/rhoB[gpI];
        const scalar k = 0.5*magSqr(U);

        heB[gpI] = EB[gpI]/rhoB[gpI] - k;

        const scalar T =
            clampTemperature(calcRawTemperatureFromRhoHe(rhoB[gpI], heB[gpI]), TMean);

        TB[gpI] = T;
        pB[gpI] = eos().calcPFromRhoT(rhoB[gpI], T);

        const scalar Cp = thermo().calcCp(T);
        const scalar Cv = thermo().calcCv(T);
        const scalar gamma = thermo().calcGamma(Cp, Cv);
        aB[gpI] = thermo().calcSpeedOfSound(T, gamma);
        muB[gpI] = transport().calcMu(T);
    }

    TG.faceField().assignBCGhostState(TB);
    pG.faceField().assignBCGhostState(pB);
    heG.faceField().assignBCGhostState(heB);
    aG.faceField().assignBCGhostState(aB);
    muG.faceField().assignBCGhostState(muB);
}

} // End namespace Foam

// ************************************************************************* //
