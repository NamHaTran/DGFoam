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
    // Build the momentum flux explicitly at each Gauss point so we can
    // benchmark the cost of GaussField operator overloading itself.
    tmp<GaussField<tensor>> tMomentumFlux = GaussField<tensor>::New(cellI, &dgMesh);
    GaussField<tensor>& momentumFlux = tMomentumFlux.ref();

    {
        const cellGaussField<vector>& rhoUCell = rhoUG.cellField();
        const cellGaussField<vector>& UCell = UG.cellField();
        const cellGaussField<scalar>& pCell = pG.cellField();
        cellGaussField<tensor>& momentumCell = momentumFlux.cellField();

        for (label gpI = 0; gpI < momentumCell.size(); ++gpI)
        {
            momentumCell[gpI] = rhoUCell[gpI]*UCell[gpI] + pCell[gpI]*symmTensor::I;
        }

        const faceGaussField<vector>& rhoUFace = rhoUG.faceField();
        const faceGaussField<vector>& UFace = UG.faceField();
        const faceGaussField<scalar>& pFace = pG.faceField();
        faceGaussField<tensor>& momentumFace = momentumFlux.faceField();

        for (label gpI = 0; gpI < momentumFace.nGauss(); ++gpI)
        {
            momentumFace.minusValueAt(gpI) =
                rhoUFace.minusValue(gpI)*UFace.minusValue(gpI)
              + pFace.minusValue(gpI)*symmTensor::I;

            momentumFace.plusValueAt(gpI) =
                rhoUFace.plusValue(gpI)*UFace.plusValue(gpI)
              + pFace.plusValue(gpI)*symmTensor::I;
        }
    }

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
    tmp<GaussField<vector>> tEnergyFlux = GaussField<vector>::New(cellI, &dgMesh);
    GaussField<vector>& energyFlux = tEnergyFlux.ref();

    {
        const cellGaussField<scalar>& ECell = EG.cellField();
        const cellGaussField<scalar>& pCell = pG.cellField();
        const cellGaussField<vector>& UCell = UG.cellField();
        cellGaussField<vector>& energyCell = energyFlux.cellField();

        for (label gpI = 0; gpI < energyCell.size(); ++gpI)
        {
            energyCell[gpI] = (ECell[gpI] + pCell[gpI])*UCell[gpI];
        }

        const faceGaussField<scalar>& EFace = EG.faceField();
        const faceGaussField<scalar>& pFace = pG.faceField();
        const faceGaussField<vector>& UFace = UG.faceField();
        faceGaussField<vector>& energyFace = energyFlux.faceField();

        for (label gpI = 0; gpI < energyFace.nGauss(); ++gpI)
        {
            energyFace.minusValueAt(gpI) =
                (EFace.minusValue(gpI) + pFace.minusValue(gpI))
               *UFace.minusValue(gpI);

            energyFace.plusValueAt(gpI) =
                (EFace.plusValue(gpI) + pFace.plusValue(gpI))
               *UFace.plusValue(gpI);
        }
    }

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
