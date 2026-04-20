//======================================================================//
//                    Total-energy conservation
//======================================================================//
//
//     d(E)/dt + div(Fc_E + Fd_E) = S_E
//
// with E storing rho*E_specific in the code and
//
//     Fc_E = (E + p) U,
//     Fd_E = tauT.U - (E/rho)mD - kappa grad(T),
//     S_E  = p tr(gradU) + tauT:gradU
//          - mD.grad(E/rho) + rho g.U.
//
// The term (E/rho) div(mD) is written as
//
//     div((E/rho)mD) - mD.grad(E/rho),
//
// so no second derivative of mD is required.
//
// Assemble the convective/pressure energy flux
//
//     Fc_E = (E + p) U.
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

// Keep the NSF molecular/thermal energy flux separate from Franz-Durst
// self-diffusion so each can select a distinct diffusive numerical flux.
const auto energyDiffFluxExpr =
    (dg::expr(tMomentumDiffFlux()) & dg::expr(UG))
  - dg::expr(kappaG)*dg::expr(gradTG);

tmp<GaussField<vector>> tEnergyDiffFlux =
    dg::evaluate<vector>(energyDiffFluxExpr, cellI, &dgMesh);

tmp<GaussField<vector>> tSelfDiffusionEnergyFlux;

if (selfDiffusion)
{
    const auto massDiffFluxExpr =
       -dg::expr(muG)
       *(
            dg::expr(gradPG)/dg::expr(pG)
          - scalar(0.5)*dg::expr(gradTG)/dg::expr(TG)
        );

    tmp<GaussField<vector>> tMassDiffFlux =
        dg::evaluate<vector>(massDiffFluxExpr, cellI, &dgMesh);

    compressibleBC->correctSelfDiffusionFlux(tMassDiffFlux.ref());

    // Assemble the Franz-Durst energy flux
    //
    //     Fd_E,self = tau_mD.U - (E/rho)mD.
    //
    // The NSF part, (-tau).U - kappa grad(T), stays in energyDiffTerm.
    const auto selfDiffusionEnergyFluxExpr =
        (dg::expr(tSelfDiffusionMomentumFlux()) & dg::expr(UG))
      - (dg::expr(EG)/dg::expr(rhoG))*dg::expr(tMassDiffFlux());

    tSelfDiffusionEnergyFlux =
        dg::evaluate<vector>
        (
            selfDiffusionEnergyFluxExpr,
            cellI,
            &dgMesh
        );

    // Assemble the gradient of the specific total energy
    //
    //     grad(E/rho) = (rho grad(E) - E grad(rho))/rho^2.
    //
    // In code notation, E is already rho*E_specific, SE = grad(E), and
    // SRho = grad(rho).
    const auto gradSpecificEnergyExpr =
        (
            dg::expr(SEG)*dg::expr(rhoG)
          - dg::expr(EG)*dg::expr(SRhoG)
        )/(dg::expr(rhoG)*dg::expr(rhoG));

    // Assemble the energy volume source
    //
    //     S_E = p tr(gradU)
    //         + tauT:gradU
    //         - mD.grad(E/rho)
    //         + rho g.U.
    const auto totalMomentumDiffFluxExpr =
        dg::expr(tMomentumDiffFlux())
      + dg::expr(tSelfDiffusionMomentumFlux());

    const auto energySourceExpr =
        dg::expr(pG)*dg::tr(dg::expr(gradUG))
      + (totalMomentumDiffFluxExpr && dg::expr(gradUG))
      - (dg::expr(tMassDiffFlux()) & gradSpecificEnergyExpr);
      // Gravity power is temporarily disabled while isolating self-diffusion:
      // + dg::expr(rhoG)*(dg::expr(UG) & dg::uniform(g.value()));

    tmp<GaussField<scalar>> tEnergySource =
        dg::evaluate<scalar>(energySourceExpr, cellI, &dgMesh);

    dgSourceTerm<scalar> energySourceTerm
    (
        "energySourceTerm",
        basisField,
        tEnergySource()
    );

    energySourceTerm.calcResidual();

    timeDiscretization.addResidual
    (
        E.name(),
        cellI,
        energySourceTerm.R()
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
        fluxSolverManager
    );

    selfDiffusionEnergyTerm.updateFlux();
    selfDiffusionEnergyTerm.calcResidual();

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
