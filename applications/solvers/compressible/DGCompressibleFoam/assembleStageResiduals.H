// Assemble the semi-discrete DG residual cell by cell.
// Each dgGeneralPDETerm computes
//     R(U) = volume contribution + surface contribution
// for one conservation law and hands the modal residual to the time integrator.
// This file is intentionally written in a "physics first" style: most flux
// expressions below are direct algebra on GaussField objects, so they read
// similarly to the continuum equations.
forAll(mesh.cells(), cellI)
{
    // dgBasisField exposes basis values/derivatives for the current cell and
    // is reused by every PDE term assembled on this cell.
    const Foam::dgBasisField& basisField = dgMesh.basisField(cellI);

    GaussField<vector>& UG      = U.gaussFields()[cellI];
    GaussField<scalar>& pG      = p.gaussFields()[cellI];
    GaussField<scalar>& EG      = E.gaussFields()[cellI];
    GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];

    // GaussField is the per-cell high-level container. It bundles:
    // - cellField(): interior quadrature values
    // - faceField(): quadrature values on faces, including ghost/neighbor data
    // dgGeneralPDETerm then consumes these Gauss-point fields and performs
    // the DG integrals for one equation on this one cell.

    // Mass equation:
    //     d(rho)/dt + div(rho U) = 0
    dgGeneralPDETerm<scalar, vector> massConvTerm
    (
        "massConvTerm",
        dgFluxSolver::equationType::massTransport,
        basisField,
        rhoG,
        rhoUG,
        fluxSolverManager
    );

    // Momentum equation:
    //     d(rhoU)/dt + div(rhoU (x) U + pI) = 0
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

    dgGeneralPDETerm<vector, tensor> momentumConvTerm
    (
        "momentumConvTerm",
        dgFluxSolver::equationType::momentumTransport,
        basisField,
        rhoUG,
        tMomentumFlux,
        fluxSolverManager
    );

    // Total-energy equation:
    //     d(E)/dt + div((E + p)U) = 0
    //
    // Again, EG + pG is only a lazy recipe here.
    // The intermediate sum is never stored as its own GaussField.
    const auto energyFluxExpr =
        (dg::expr(EG) + dg::expr(pG))*dg::expr(UG);

    tmp<GaussField<vector>> tEnergyFlux =
        dg::evaluate<vector>(energyFluxExpr, cellI, &dgMesh);

    dgGeneralPDETerm<scalar, vector> energyConvTerm
    (
        "energyConvTerm",
        dgFluxSolver::equationType::energyTransport,
        basisField,
        EG,
        tEnergyFlux,
        fluxSolverManager
    );

    // updateFlux() computes the numerical surface fluxes requested by the
    // term's configured flux solver; calcResidual() then combines volume and
    // surface contributions into the modal residual vector R().
    // In other words:
    // - dgGeneralPDETerm knows how to integrate a DG term,
    // - dgFluxSolverManager knows how to close it with a numerical flux.
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

    energyConvTerm.updateFlux();
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
        E.name(),
        cellI,
        energyConvTerm.R()
    );
}
