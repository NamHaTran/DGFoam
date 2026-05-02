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

momentumConvTerm.add(rhoUG, tMomentumFlux);

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

momentumDiffTerm.add(rhoUG, tMomentumDiffFlux);

momentumConvTerm.updateFlux(cellI);
momentumDiffTerm.updateFlux(cellI);

dgResidual<vector, tensor> momentumResidual(dgMesh, cellI);
momentumResidual.addTerm(momentumConvTerm);
momentumResidual.addTerm(momentumDiffTerm);

tmp<GaussField<tensor>> tArtificialDiffusionMomentumFlux;

if (artificialDiffusion)
{
    const auto artificialDiffusionMomentumFluxExpr =
       -dg::expr(muAVG)*dg::expr(SRhoUG);

    tArtificialDiffusionMomentumFlux =
        dg::evaluate<tensor>
        (
            artificialDiffusionMomentumFluxExpr,
            cellI,
            &dgMesh
        );

    artificialDiffusionMomentumTerm->add
    (
        rhoUG,
        tArtificialDiffusionMomentumFlux
    );

    artificialDiffusionMomentumTerm->updateFlux(cellI);
    momentumResidual.addTerm(*artificialDiffusionMomentumTerm);
}

momentumResidual.calcResidual();

timeDiscretization.addResidual
(
    rhoU.name(),
    cellI,
    momentumResidual.R()
);

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
//         dgMesh,
//         tMomentumSource()
//     );
//
//     momentumResidual.addTerm(momentumSourceTerm);
//
//     timeDiscretization.addResidual
//     (
//         rhoU.name(),
//         cellI,
//         momentumResidual.R()
//     );
// }
