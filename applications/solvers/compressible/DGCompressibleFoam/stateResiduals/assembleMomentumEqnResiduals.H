//======================================================================//
//                      Momentum conservation
//======================================================================//
//
//     d(rhoU)/dt + div(Fc_m + Fd_m) = rho*g
//
// with
//
//     Fc_m = rhoU (x) U + pI,
//     Fd_m = tauT,
//     tauT = -tau + mD (x) U + U (x) mD - 2/3 (mD.U) I,
//     tau  = mu*(gradU + gradU^T - 2/3 tr(gradU) I).
//
// Assemble the convective/pressure momentum flux
//
//     Fc_m = rhoU (x) U + pI.
const auto momentumFluxExpr =
    dg::expr(rhoUG)*dg::expr(UG)
  + dg::expr(pG)*dg::uniform(symmTensor::I);

tmp<GaussField<tensor>> tMomentumFlux =
    dg::evaluate<tensor>(momentumFluxExpr, cellI, &dgMesh);

dgConvectiveTerm<vector, tensor> momentumConvTerm
(
    "momentumConvTerm",
    dgFluxSolver::equationType::momentumTransport,
    dgFluxSolver::fluxType::convective,
    basisField,
    rhoUG,
    tMomentumFlux,
    fluxSolverManager
);

// Assemble the Newtonian viscous stress
//
//     tau = mu*(gradU + gradU^T - 2/3 tr(gradU) I).
//
// The conservative NSF diffusive flux uses -tau.  When Franz-Durst
// self-diffusion is enabled, add its molecular-momentum flux before
// constructing the single momentum diffusive term:
//
//     Fd_m = -tau + tau_mD,
//
// with
//
//     tau_mD = mD (x) U + U (x) mD - 2/3 (mD.U) I.
const auto tauExpr =
    dg::expr(muG)
   *(
        dg::expr(gradUG)
      + dg::T(dg::expr(gradUG))
      - scalar(2.0/3.0)
       *dg::tr(dg::expr(gradUG))
       *dg::uniform(symmTensor::I)
    );

tmp<GaussField<tensor>> tMomentumDiffFlux;

if (selfDiffusion)
{
    const auto massMomentumDiffFluxExpr =
        dg::expr(tMassDiffFlux())*dg::expr(UG)
      + dg::expr(UG)*dg::expr(tMassDiffFlux())
      - scalar(2.0/3.0)
       *(dg::expr(tMassDiffFlux()) & dg::expr(UG))
       *dg::uniform(symmTensor::I);

    const auto totalMomentumDiffFluxExpr =
        -tauExpr + massMomentumDiffFluxExpr;

    tMomentumDiffFlux =
        dg::evaluate<tensor>
        (
            totalMomentumDiffFluxExpr,
            cellI,
            &dgMesh
        );
}
else
{
    tMomentumDiffFlux =
        dg::evaluate<tensor>
        (
            -tauExpr,
            cellI,
            &dgMesh
        );
}

dgDiffusiveTerm<vector, tensor> momentumDiffTerm
(
    "momentumDiffTerm",
    dgFluxSolver::equationType::momentumTransport,
    dgFluxSolver::fluxType::diffusive,
    basisField,
    rhoUG,
    tMomentumDiffFlux,
    fluxSolverManager
);

momentumConvTerm.updateFlux();
momentumDiffTerm.updateFlux();
momentumConvTerm.calcResidual();
momentumDiffTerm.calcResidual();

timeDiscretization.addResidual
(
    rhoU.name(),
    cellI,
    momentumConvTerm.R()
);

timeDiscretization.addResidual
(
    rhoU.name(),
    cellI,
    momentumDiffTerm.R()
);

if (artificialDiffusion)
{
    const auto artificialDiffusionMomentumFluxExpr =
       -dg::expr(muAVG)*dg::expr(SRhoUG);

    tmp<GaussField<tensor>> tArtificialDiffusionMomentumFlux =
        dg::evaluate<tensor>
        (
            artificialDiffusionMomentumFluxExpr,
            cellI,
            &dgMesh
        );

    dgDiffusiveTerm<vector, tensor> artificialDiffusionMomentumTerm
    (
        "artificialDiffusionMomentumTerm",
        dgFluxSolver::equationType::momentumTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        rhoUG,
        tArtificialDiffusionMomentumFlux,
        fluxSolverManager
    );

    artificialDiffusionMomentumTerm.updateFlux();
    artificialDiffusionMomentumTerm.calcResidual();

    timeDiscretization.addResidual
    (
        rhoU.name(),
        cellI,
        artificialDiffusionMomentumTerm.R()
    );
}

// Gravity source is temporarily disabled for both selfDiffusion on/off while
// isolating the self-diffusion terms.
//
// if (selfDiffusion)
// {
//     //     S_m = rho*g.
//     const auto momentumSourceExpr =
//         dg::expr(rhoG)*dg::uniform(g.value());
//
//     tmp<GaussField<vector>> tMomentumSource =
//         dg::evaluate<vector>(momentumSourceExpr, cellI, &dgMesh);
//
//     dgSourceTerm<vector> momentumSourceTerm
//     (
//         "momentumSourceTerm",
//         basisField,
//         tMomentumSource()
//     );
//
//     momentumSourceTerm.calcResidual();
//
//     timeDiscretization.addResidual
//     (
//         rhoU.name(),
//         cellI,
//         momentumSourceTerm.R()
//     );
// }
