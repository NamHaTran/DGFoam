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
    dgFluxSolver::fluxType::convective,
    basisField,
    rhoG,
    rhoUG,
    fluxSolverManager
);

massConvTerm.updateFlux();
massConvTerm.calcResidual();

if (selfDiffusion)
{
    // Franz-Durst self-diffusive mass flux mD is preassembled once per cell
    // in assembleStageResiduals.H and reused across mass/momentum/energy.

    dgDiffusiveTerm<scalar, vector> selfDiffusionMassTerm
    (
        "selfDiffusionMassTerm",
        dgFluxSolver::equationType::massTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        rhoG,
        tMassDiffFlux,
        fluxSolverManager
    );

    selfDiffusionMassTerm.updateFlux();
    selfDiffusionMassTerm.calcResidual();

    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        massConvTerm.R()
    );

    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        selfDiffusionMassTerm.R()
    );
}
else
{
    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        massConvTerm.R()
    );
}

if (artificialDiffusion)
{
    const auto artificialDiffusionMassFluxExpr =
       -dg::expr(muAVG)*dg::expr(SRhoG);

    tmp<GaussField<vector>> tArtificialDiffusionMassFlux =
        dg::evaluate<vector>
        (
            artificialDiffusionMassFluxExpr,
            cellI,
            &dgMesh
        );

    dgDiffusiveTerm<scalar, vector> artificialDiffusionMassTerm
    (
        "artificialDiffusionMassTerm",
        dgFluxSolver::equationType::massTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        rhoG,
        tArtificialDiffusionMassFlux,
        fluxSolverManager
    );

    artificialDiffusionMassTerm.updateFlux();
    artificialDiffusionMassTerm.calcResidual();

    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        artificialDiffusionMassTerm.R()
    );
}
