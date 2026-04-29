//======================================================================//
//                    Total-energy conservation
//======================================================================//
//
//     d(rho Ec)/dt + div(rho U Ec)
//       = div(qT) + Ek div(mD) - div(P U) - div(tauT.U)
//
// with E storing rho*Ec in the code and
//
//     tauT = -mu*(gradU + gradU^T - 2/3 tr(gradU) I)
//          + mD (x) U + U (x) mD - 2/3 (mD.U) I,
//     qT   = -kappa grad(T) + mD e,
//     Ec   = e + Ek,
//     Ek   = 0.5*(U.U).
//
// In conservative left-hand-side form the divergence terms are assembled as
//
//     d(E)/dt + div(E U + P U + tauT.U + qT) = Ek div(mD).
//
// The kinetic mass-diffusion source is first moved to the left-hand side:
//
//     d(E)/dt + ... - Ek div(mD) = 0.
//
// Then use the product rule
//
//     Ek div(mD) = div(Ek mD) - mD.grad(Ek),
//
// giving the left-hand-side contribution
//
//     -Ek div(mD) = -div(Ek mD) + mD.grad(Ek).
//
// so for each test function phi on a cell K,
//
//     int_K phi [-Ek div(mD)] dV
//       = -int_K phi div(Ek mD) dV
//        +int_K phi mD.grad(Ek) dV
//
//       = -int_dK phi (Ek mD)^*.n dS
//        +int_K Ek mD.grad(phi) dV
//        +int_K phi mD.grad(Ek) dV.
//
// The dgDiffusiveTerm convention returns
//
//     R(F) = int_K F.grad(phi) dV - int_dK phi F^*.n dS,
//
// so choosing F = -Ek*mD assembles the -div(Ek*mD) part.  The remaining
// volume correction is added separately as +mD.grad(Ek).
//
// Assemble the convective energy flux
//
//     Fc_E = E U + P U.
const auto energyFluxExpr =
    (dg::expr(EG) + dg::expr(pG))*dg::expr(UG);

tmp<GaussField<vector>> tEnergyFlux =
    dg::evaluate<vector>(energyFluxExpr, cellI, &dgMesh);

dgConvectiveTerm<scalar, vector> energyConvTerm
(
    "energyConvTerm",
    dgFluxSolver::equationType::energyTransport,
    dgFluxSolver::fluxType::convective,
    basisField,
    EG,
    tEnergyFlux,
    fluxSolverManager
);

// Keep the pressure/NSF molecular/thermal energy flux separate from
// Franz-Durst self-diffusion so each can select a distinct diffusive
// numerical flux.
//
//     F_E,NSF = tau_NSF.U + q_NSF
//             = tau_NSF.U - kappa grad(T)
//
// where tau_NSF is stored in tMomentumDiffFlux.
const auto energyDiffFluxExpr =
  (dg::expr(tMomentumDiffFlux()) & dg::expr(UG))
  - dg::expr(kappaG)*dg::expr(gradTG);

tmp<GaussField<vector>> tEnergyDiffFlux =
    dg::evaluate<vector>(energyDiffFluxExpr, cellI, &dgMesh);

tmp<GaussField<vector>> tSelfDiffusionEnergyFlux;

if (selfDiffusion)
{
    // Assemble the Franz-Durst molecular-momentum and heat-flux part:
    //
    //     Fd_E,self = tau_mD.U + e*mD.
    //
    // Use the thermo model for e(T) instead of reconstructing it from
    // Ec - Ek.  This keeps the heat flux consistent with the active
    // caloric relation.
    tmp<GaussField<scalar>> tInternalEnergy =
        GaussField<scalar>::New(cellI, &dgMesh);

    thermo->energyModel().calcHe(cellI, TG, tInternalEnergy.ref());

    const auto selfDiffusionEnergyFluxExpr =
        (dg::expr(tSelfDiffusionMomentumFlux()) & dg::expr(UG))
      + dg::expr(tInternalEnergy())*dg::expr(tMassDiffFlux());

    tSelfDiffusionEnergyFlux =
        dg::evaluate<vector>
        (
            selfDiffusionEnergyFluxExpr,
            cellI,
            &dgMesh
        );
}

dgDiffusiveTerm<scalar, vector> energyDiffTerm
(
    "energyDiffTerm",
    dgFluxSolver::equationType::energyTransport,
    dgFluxSolver::fluxType::diffusive,
    basisField,
    EG,
    tEnergyDiffFlux,
    fluxSolverManager
);

energyConvTerm.updateFlux();
energyDiffTerm.updateFlux();
energyConvTerm.calcResidual();
energyDiffTerm.calcResidual();

if (selfDiffusion)
{
    dgDiffusiveTerm<scalar, vector> selfDiffusionEnergyTerm
    (
        "selfDiffusionEnergyTerm",
        dgFluxSolver::equationType::energyTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        EG,
        tSelfDiffusionEnergyFlux,
        fluxSolverManager,
        false
    );

    selfDiffusionEnergyTerm.updateFlux();
    selfDiffusionEnergyTerm.calcResidual();

    const auto kineticEnergyExpr =
        scalar(0.5)*(dg::expr(UG) & dg::expr(UG));

    const auto kineticMassDiffFluxExpr =
       -kineticEnergyExpr*dg::expr(tMassDiffFlux());

    tmp<GaussField<vector>> tKineticMassDiffFlux =
        dg::evaluate<vector>
        (
            kineticMassDiffFluxExpr,
            cellI,
            &dgMesh
        );

    dgDiffusiveTerm<scalar, vector> kineticMassDiffFluxTerm
    (
        "kineticMassDiffFluxTerm",
        dgFluxSolver::equationType::energyTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        EG,
        tKineticMassDiffFlux,
        fluxSolverManager,
        false
    );

    kineticMassDiffFluxTerm.updateFlux();
    kineticMassDiffFluxTerm.calcResidual();

    const auto gradKineticEnergyExpr =
        dg::expr(UG) & dg::expr(gradUG);

    const auto kineticMassDiffCorrectionExpr =
       (dg::expr(tMassDiffFlux()) & gradKineticEnergyExpr);

    tmp<GaussField<scalar>> tKineticMassDiffCorrection =
        dg::evaluate<scalar>
        (
            kineticMassDiffCorrectionExpr,
            cellI,
            &dgMesh
        );

    dgSourceTerm<scalar> kineticMassDiffCorrectionTerm
    (
        "kineticMassDiffCorrectionTerm",
        basisField,
        tKineticMassDiffCorrection(),
        false
    );

    kineticMassDiffCorrectionTerm.calcResidual();

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        energyConvTerm.R()
    );

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        energyDiffTerm.R()
    );

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        selfDiffusionEnergyTerm.R()
    );

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        kineticMassDiffFluxTerm.R()
    );

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        kineticMassDiffCorrectionTerm.R()
    );
}
else
{
    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        energyConvTerm.R()
    );

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        energyDiffTerm.R()
    );
}

if (artificialDiffusion)
{
    const auto artificialDiffusionEnergyFluxExpr =
       -dg::expr(muAVG)*dg::expr(SEG);

    tmp<GaussField<vector>> tArtificialDiffusionEnergyFlux =
        dg::evaluate<vector>
        (
            artificialDiffusionEnergyFluxExpr,
            cellI,
            &dgMesh
        );

    dgDiffusiveTerm<scalar, vector> artificialDiffusionEnergyTerm
    (
        "artificialDiffusionEnergyTerm",
        dgFluxSolver::equationType::energyTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        EG,
        tArtificialDiffusionEnergyFlux,
        fluxSolverManager
    );

    artificialDiffusionEnergyTerm.updateFlux();
    artificialDiffusionEnergyTerm.calcResidual();

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        artificialDiffusionEnergyTerm.R()
    );
}
