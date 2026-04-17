// Assemble the semi-discrete compressible Navier-Stokes-Fourier residual
// cell by cell.
//
// The conservative system is written in divergence form as
//
//     dQ/dt + div(Fc(Q) + Fd(Q, gradQ)) = 0,
//
// with Q = (rho, rhoU, E).  This readable version builds the fluxes using
// dgExpr algebra on GaussField objects, so each block keeps the same shape as
// the continuum equation:
//
//     mass     : convective only
//     momentum : convective + viscous stress
//     energy   : convective + viscous work + Fourier heat conduction
forAll(mesh.cells(), cellI)
{
    // dgBasisField exposes basis values/derivatives for the current cell and
    // is reused by every PDE term assembled on this cell.
    const Foam::dgBasisField& basisField = dgMesh.basisField(cellI);

    GaussField<vector>& UG      = U.gaussFields()[cellI];
    GaussField<scalar>& TG      = T.gaussFields()[cellI];
    GaussField<scalar>& pG      = p.gaussFields()[cellI];
    GaussField<scalar>& EG      = E.gaussFields()[cellI];
    GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];
    GaussField<tensor>& gradUG  = gradU.gaussFields()[cellI];
    GaussField<vector>& gradTG  = gradT.gaussFields()[cellI];

    // GaussField is the per-cell high-level container. It bundles:
    // - cellField(): interior quadrature values
    // - faceField(): quadrature values on faces, including ghost/neighbor data
    // dgConvectiveTerm then consumes these Gauss-point fields and performs
    // the DG integrals for one equation on this one cell.

    //==================================================================//
    //                          Mass conservation
    //==================================================================//
    //
    // General equation:
    //
    //     d(rho)/dt + div(Fc_rho) = 0
    //
    // with
    //
    //     Fc_rho = rho U.
    //
    // There is no Navier-Stokes-Fourier diffusive mass flux in this model.
    //
    // Convective term assembly:
    // rhoUG already stores rho*U at cell and face Gauss points.
    dgConvectiveTerm<scalar, vector> massConvTerm
    (
        "massConvTerm",
        dgFluxSolver::equationType::massTransport,
        basisField,
        rhoG,
        rhoUG,
        fluxSolverManager
    );

    //==================================================================//
    //                       Momentum conservation
    //==================================================================//
    //
    // General equation:
    //
    //     d(rhoU)/dt + div(Fc_m + Fd_m) = 0
    //
    // with
    //
    //     Fc_m = rhoU (x) U + pI,
    //     Fd_m = -tau,
    //     tau  = mu*(gradU + gradU^T - 2/3 tr(gradU) I).
    //
    // Momentum convective term assembly:
    // Build Fc_m at all Gauss points, then let dgConvectiveTerm integrate it
    // with the fluxSolverManager-selected numerical flux.
    //
    // This is the minimal dgExpr prototype in action.
    //
    // Important reading guide:
    // 1. dg::expr(field) does NOT compute anything.
    //    It just creates a tiny leaf node that references the field.
    // 2. operator* and operator+ below do NOT return a tmp<GaussField>.
    //    They create a lightweight expression tree.
    // 3. dg::evaluate<tensor>(...) is the point where the tree is finally
    //    walked and materialized into a real GaussField storage object.
    //
    // So this line keeps the mathematical form easy to read, but the actual
    // loop over Gauss points is hidden inside dg::assign()/dg::evaluate().
    const auto momentumFluxExpr =
        dg::expr(rhoUG)*dg::expr(UG)
      + dg::expr(pG)*dg::uniform(symmTensor::I);

    tmp<GaussField<tensor>> tMomentumFlux =
        dg::evaluate<tensor>(momentumFluxExpr, cellI, &dgMesh);

    dgConvectiveTerm<vector, tensor> momentumConvTerm
    (
        "momentumConvTerm",
        dgFluxSolver::equationType::momentumTransport,
        basisField,
        rhoUG,
        tMomentumFlux,
        fluxSolverManager
    );

    // The transport properties are evaluated through the selected thermo
    // transport model, so Sutherland/power-law/etc. stay runtime-selected.
    // muG and kappaG are full GaussField objects because the diffusive flux
    // needs both cell values and face minus/plus values.
    GaussField<scalar> muG(cellI, &dgMesh);
    GaussField<scalar> kappaG(cellI, &dgMesh);

    thermo->transport().calcMu(cellI, TG, muG);
    thermo->transport().calcKappa(cellI, TG, kappaG);

    // Momentum diffusive term assembly:
    // Build the Newtonian stress and store the conservative diffusive flux
    // Fd_m = -tau.  dgDiffusiveTerm applies the central flux on faces.
    const auto tauExpr =
        dg::expr(muG)
       *(
            dg::expr(gradUG)
          + dg::T(dg::expr(gradUG))
          - scalar(2.0/3.0)
           *dg::tr(dg::expr(gradUG))
           *dg::uniform(symmTensor::I)
        );

    tmp<GaussField<tensor>> tMomentumDiffFlux =
        dg::evaluate<tensor>(-tauExpr, cellI, &dgMesh);

    dgDiffusiveTerm<vector, tensor> momentumDiffTerm
    (
        "momentumDiffTerm",
        dgFluxSolver::equationType::momentumTransport,
        basisField,
        rhoUG,
        tMomentumDiffFlux
    );

    //==================================================================//
    //                     Total-energy conservation
    //==================================================================//
    //
    // General equation:
    //
    //     d(E)/dt + div(Fc_E + Fd_E) = 0
    //
    // with
    //
    //     Fc_E = (E + p) U,
    //     Fd_E = -(tau.U + kappa grad(T)).
    //
    // Energy convective term assembly:
    // EG + pG is only a lazy expression recipe here; the intermediate sum is
    // never stored as its own GaussField.
    const auto energyFluxExpr =
        (dg::expr(EG) + dg::expr(pG))*dg::expr(UG);

    tmp<GaussField<vector>> tEnergyFlux =
        dg::evaluate<vector>(energyFluxExpr, cellI, &dgMesh);

    dgConvectiveTerm<scalar, vector> energyConvTerm
    (
        "energyConvTerm",
        dgFluxSolver::equationType::energyTransport,
        basisField,
        EG,
        tEnergyFlux,
        fluxSolverManager
    );

    // Energy diffusive term assembly:
    // The first part is viscous work, the second is Fourier conduction.
    const auto energyDiffFluxExpr =
       -(
            (tauExpr & dg::expr(UG))
          + dg::expr(kappaG)*dg::expr(gradTG)
        );

    tmp<GaussField<vector>> tEnergyDiffFlux =
        dg::evaluate<vector>(energyDiffFluxExpr, cellI, &dgMesh);

    dgDiffusiveTerm<scalar, vector> energyDiffTerm
    (
        "energyDiffTerm",
        dgFluxSolver::equationType::energyTransport,
        basisField,
        EG,
        tEnergyDiffFlux
    );

    //==================================================================//
    //                         Residual assembly
    //==================================================================//
    //
    // updateFlux() computes the face flux:
    // - convective terms use fluxSolverManager-selected numerical fluxes,
    // - diffusive terms use dgDiffusiveTerm's central flux.
    //
    // calcResidual() then combines volume and surface contributions into the
    // modal residual vector R() for each conservation law.
    massConvTerm.updateFlux();
    // debugProcessorFaceExchange::printFaceFlux
    // (
    //     "mass numerical flux",
    //     massConvTerm.fluxField()
    // );
    // debugProcessorFaceExchange::printFaceIntegral
    // (
    //     "mass face integral",
    //     massConvTerm.fluxField()
    // );
    massConvTerm.calcResidual();

    momentumConvTerm.updateFlux();
    momentumDiffTerm.updateFlux();
    // debugProcessorFaceExchange::printFaceFlux
    // (
    //     "momentum numerical flux",
    //     momentumConvTerm.fluxField()
    // );
    // debugProcessorFaceExchange::printFaceIntegral
    // (
    //     "momentum face integral",
    //     momentumConvTerm.fluxField()
    // );
    momentumConvTerm.calcResidual();
    momentumDiffTerm.calcResidual();

    energyConvTerm.updateFlux();
    energyDiffTerm.updateFlux();
    // debugProcessorFaceExchange::printFaceFlux
    // (
    //     "energy numerical flux",
    //     energyConvTerm.fluxField()
    // );
    // debugProcessorFaceExchange::printFaceIntegral
    // (
    //     "energy face integral",
    //     energyConvTerm.fluxField()
    // );
    energyConvTerm.calcResidual();
    energyDiffTerm.calcResidual();

    //==================================================================//
    //                   Store residuals for time update
    //==================================================================//
    //
    // The time discretization object stores residuals per registered field
    // and per stage. The actual update is deferred to advanceStage.H.
    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        massConvTerm.R()
    );

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
