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
massConvTerm.add(rhoG, rhoUG);
massConvTerm.updateFlux(cellI);

dgResidual<scalar, vector> massResidual(dgMesh, cellI);
massResidual.addTerm(massConvTerm);

if (selfDiffusion)
{
    // Franz-Durst self-diffusive mass flux mD is preassembled once per cell
    // in assembleStageResiduals.H and reused across mass/momentum/energy.
    selfDiffusionMassTerm->add(rhoG, tMassDiffFlux);
    selfDiffusionMassTerm->updateFlux(cellI);
    massResidual.addTerm(*selfDiffusionMassTerm);
}

tmp<GaussField<vector>> tArtificialDiffusionMassFlux;

if (artificialDiffusion)
{
    const auto artificialDiffusionMassFluxExpr =
       -dg::expr(muAVG)*dg::expr(SRhoG);

    tArtificialDiffusionMassFlux =
        dg::evaluate<vector>
        (
            artificialDiffusionMassFluxExpr,
            cellI,
            &dgMesh
        );

    artificialDiffusionMassTerm->add(rhoG, tArtificialDiffusionMassFlux);
    artificialDiffusionMassTerm->updateFlux(cellI);
    massResidual.addTerm(*artificialDiffusionMassTerm);
}

massResidual.calcResidual();

timeDiscretization.addResidual
(
    rho.name(),
    cellI,
    massResidual.R()
);
