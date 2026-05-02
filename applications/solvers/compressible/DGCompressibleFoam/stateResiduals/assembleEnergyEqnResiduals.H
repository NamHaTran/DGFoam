//======================================================================//
//                    Total-energy conservation
//======================================================================//
//
//     d(rho Ec)/dt + div(rho U Ec)
//       = div(qT) + Ek div(mD) - div(P U) - div(tauT.U)
//         + rho U.g
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
// Gravity contributes the total-energy source
//
//     S_E = rho U.g = rhoU.g.
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

// Assemble the pressure/NSF molecular/thermal energy flux.  Since
// tMomentumDiffFlux already contains the total momentum diffusive flux,
// including Franz-Durst self-diffusion when enabled, this term contains
// tauT.U directly:
//
//     F_E = tauT.U - kappa grad(T)
//
// When self-diffusion is active, add the heat/mass and conservative kinetic
// mass-diffusion fluxes before constructing the single energy diffusive term:
//
//     F_E,self = (e - Ek)*mD.
tmp<GaussField<vector>> tEnergyDiffFlux;

if (selfDiffusion)
{
    // Use the thermo model for e(T) instead of reconstructing it from
    // Ec - Ek.  This keeps the heat flux consistent with the active
    // caloric relation.
    tmp<GaussField<scalar>> tInternalEnergy =
        GaussField<scalar>::New(cellI, &dgMesh);

    thermo->energyModel().calcHe(cellI, TG, tInternalEnergy.ref());

    const auto kineticEnergyExpr =
        scalar(0.5)*(dg::expr(UG) & dg::expr(UG));

    const auto totalEnergyDiffFluxExpr =
        (dg::expr(tMomentumDiffFlux()) & dg::expr(UG))
      - dg::expr(kappaG)*dg::expr(gradTG)
      + (dg::expr(tInternalEnergy()) - kineticEnergyExpr)
       *dg::expr(tMassDiffFlux());

    tEnergyDiffFlux =
        dg::evaluate<vector>
        (
            totalEnergyDiffFluxExpr,
            cellI,
            &dgMesh
        );
}
else
{
    const auto energyDiffFluxExpr =
        (dg::expr(tMomentumDiffFlux()) & dg::expr(UG))
      - dg::expr(kappaG)*dg::expr(gradTG);

    tEnergyDiffFlux =
        dg::evaluate<vector>(energyDiffFluxExpr, cellI, &dgMesh);
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

dgResidual<scalar, vector> energyResidual(basisField);
energyResidual.addTerm(energyConvTerm);
energyResidual.addTerm(energyDiffTerm);

std::optional<dgSourceTerm<scalar>> kineticMassDiffCorrectionTerm;

if (selfDiffusion)
{
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

    kineticMassDiffCorrectionTerm.emplace
    (
        "kineticMassDiffCorrectionTerm",
        basisField,
        tKineticMassDiffCorrection()
    );

    energyResidual.addTerm(*kineticMassDiffCorrectionTerm);
}

std::optional<dgDiffusiveTerm<scalar, vector>>
    artificialDiffusionEnergyTerm;

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

    artificialDiffusionEnergyTerm.emplace
    (
        "artificialDiffusionEnergyTerm",
        dgFluxSolver::equationType::energyTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        EG,
        tArtificialDiffusionEnergyFlux,
        fluxSolverManager
    );

    artificialDiffusionEnergyTerm->updateFlux();
    energyResidual.addTerm(*artificialDiffusionEnergyTerm);
}

//     S_E = rhoU.g.
const auto energyGravitySourceExpr =
    dg::expr(rhoUG) & dg::uniform(g.value());

tmp<GaussField<scalar>> tEnergyGravitySource =
    dg::evaluate<scalar>(energyGravitySourceExpr, cellI, &dgMesh);

dgSourceTerm<scalar> energyGravitySourceTerm
(
    "energyGravitySourceTerm",
    basisField,
    tEnergyGravitySource()
);

energyResidual.addTerm(energyGravitySourceTerm);

energyResidual.calcResidual();

timeDiscretization.addResidual
(
    E.name(),
    cellI,
    energyResidual.R()
);
