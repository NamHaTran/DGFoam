// Assemble the semi-discrete DG residuals for all conservative equations.
for (label cellI = 0; cellI < mesh.C().size(); cellI++)
{
    Foam::dgBasisField basisField(cellI, dgMesh);

    GaussField<vector>& UG      = U.gaussFields()[cellI];
    GaussField<scalar>& pG      = p.gaussFields()[cellI];
    GaussField<scalar>& EG      = E.gaussFields()[cellI];
    GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];

    dgGeneralPDETerm<scalar, vector> massConvTerm
    (
        "massConvTerm",
        dgFluxSolver::equationType::massTransport,
        basisField,
        rhoG,
        rhoUG,
        fluxSolverManager
    );

    // Inviscid momentum equation:
    // d(rhoU)/dt + div(rhoU (x) U + pI) = 0
    GaussField<symmTensor> identityTensor
    (
        cellI,
        &dgMesh,
        symmTensor::I
    );

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

    // Inviscid total-energy equation:
    // d(E)/dt + div((E + p)U) = 0
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

    massConvTerm.updateFlux();
    massConvTerm.calcResidual();

    momentumConvTerm.updateFlux();
    momentumConvTerm.calcResidual();

    energyConvTerm.updateFlux();
    energyConvTerm.calcResidual();

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

    Pout<< "Stage " << timeDiscretization.stage()
        << " cell " << cellI
        << " mass residual: " << massConvTerm.R()
        << " momentum residual: " << momentumConvTerm.R()
        << " energy residual: " << energyConvTerm.R() << endl;
}
