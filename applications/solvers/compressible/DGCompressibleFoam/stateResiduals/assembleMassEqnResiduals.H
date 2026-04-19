//======================================================================//
//                         Mass conservation
//======================================================================//
//
//     d(rho)/dt + div(Fc_rho + Fd_rho) = 0
//
// with
//
//     Fc_rho = rho U,
//     Fd_rho = mD
//            = -mu*(grad(p)/p - 0.5*grad(T)/T).
//
// The pressure and temperature gradients were reconstructed by
// solveAuxiliaryEqn.H before this residual assembly.
dgConvectiveTerm<scalar, vector> massConvTerm
(
    "massConvTerm",
    dgFluxSolver::equationType::massTransport,
    basisField,
    rhoG,
    rhoUG,
    fluxSolverManager
);

massConvTerm.updateFlux();
massConvTerm.calcResidual();

timeDiscretization.addResidual
(
    rho.name(),
    cellI,
    massConvTerm.R()
);

if (selfDiffusion)
{
    // Assemble the Franz-Durst self-diffusive mass flux
    //
    //     Fd_rho = mD
    //            = -mu*(grad(p)/p - 0.5*grad(T)/T).
    //
    // This flux enters the mass equation as
    //
    //     d(rho)/dt + div(rho U + mD) = 0.
    const auto massDiffFluxExpr =
       -dg::expr(muG)
       *(
            dg::expr(gradPG)/dg::expr(pG)
          - scalar(0.5)*dg::expr(gradTG)/dg::expr(TG)
        );

    tmp<GaussField<vector>> tMassDiffFlux =
        dg::evaluate<vector>(massDiffFluxExpr, cellI, &dgMesh);

    compressibleBC->correctSelfDiffusionFlux(tMassDiffFlux.ref());

    dgDiffusiveTerm<scalar, vector> massDiffTerm
    (
        "massDiffTerm",
        dgFluxSolver::equationType::massTransport,
        basisField,
        rhoG,
        tMassDiffFlux
    );

    massDiffTerm.updateFlux();
    massDiffTerm.calcResidual();

    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        massDiffTerm.R()
    );
}
