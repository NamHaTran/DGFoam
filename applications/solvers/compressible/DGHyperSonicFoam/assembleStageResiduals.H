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
    // The identity tensor is promoted to a Gauss field so pI can be formed
    // pointwise at all quadrature points.
    GaussField<symmTensor> identityTensor
    (
        cellI,
        &dgMesh,
        symmTensor::I
    );

    // Thanks to operator overloading, rhoUG*UG and pG*identityTensor are
    // pointwise Gauss-field products. This is one of the main reasons DG
    // equation assembly in dgFoam stays close to the mathematical notation.
    tmp<GaussField<tensor>> tMomentumFlux = rhoUG*UG + pG*identityTensor;

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
    tmp<GaussField<vector>> tEnergyFlux = (EG + pG)*UG;

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
    massConvTerm.calcResidual();

    momentumConvTerm.updateFlux();
    momentumConvTerm.calcResidual();

    energyConvTerm.updateFlux();
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
