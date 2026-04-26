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
// The div(mD) operator is handled by a dgDiffusiveTerm, as in the mass
// equation.
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
        true
    );

    selfDiffusionEnergyTerm.updateFlux();
    selfDiffusionEnergyTerm.calcResidual();

    dgDiffusiveTerm<scalar, vector> selfDiffusionMassDivergenceTerm
    (
        "selfDiffusionMassTerm",
        dgFluxSolver::equationType::massTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        rhoG,
        tMassDiffFlux,
        fluxSolverManager,
        true
    );

    selfDiffusionMassDivergenceTerm.updateFlux();
    selfDiffusionMassDivergenceTerm.calcResidual();

    List<scalar> kineticMassDiffResidual
    (
        selfDiffusionMassDivergenceTerm.R().size(),
        pTraits<scalar>::zero
    );

    const cellDof<scalar>& rhoDof = rho.dof()[cellI];
    const cellDof<vector>& rhoUDof = rhoU.dof()[cellI];

    forAll(kineticMassDiffResidual, dofI)
    {
        const scalar rhoValue =
            max(mag(rhoDof[dofI]), VSMALL);
        const vector UValue = rhoUDof[dofI]/rhoValue;

        const scalar kineticEnergy =
            scalar(0.5)*(UValue & UValue);

        // selfDiffusionMassDivergenceTerm.R() is the residual contribution
        // for +div(mD) on the left-hand side of the mass equation, i.e.
        // -div(mD) in the explicit update.  The energy equation has
        // +Ek*div(mD) on the right-hand side, so the sign is reversed here.
        kineticMassDiffResidual[dofI] =
           -kineticEnergy*selfDiffusionMassDivergenceTerm.R()[dofI];
    }

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
        kineticMassDiffResidual
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
