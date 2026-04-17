// Optimized semi-discrete compressible Navier-Stokes-Fourier residual
// assembly.
//
// This version mirrors assembleStageResiduals.H but builds algebraic flux
// fields with explicit loops over Gauss points. That is a little more verbose,
// but it can improve performance because it avoids constructing temporary
// expression/evaluation fields for these pointwise fluxes.
//
// The conservative system is assembled equation by equation:
//
//     mass     : d(rho)/dt  + div(Fc_rho)        = 0
//     momentum : d(rhoU)/dt + div(Fc_m + Fd_m)   = 0
//     energy   : d(E)/dt    + div(Fc_E + Fd_E)   = 0
//
// where the Navier-Stokes-Fourier diffusive fluxes are
//
//     Fd_m = -tau,
//     Fd_E = -(tau.U + kappa grad(T)).
//
// Use assembleStageResiduals.H when readability and compact mathematical
// expressions are preferred. Use this optimized include when the cost of
// temporary field construction matters more.
forAll(mesh.cells(), cellI)
{
    const Foam::dgBasisField& basisField = dgMesh.basisField(cellI);

    GaussField<vector>& UG      = U.gaussFields()[cellI];
    GaussField<scalar>& TG      = T.gaussFields()[cellI];
    GaussField<scalar>& pG      = p.gaussFields()[cellI];
    GaussField<scalar>& EG      = E.gaussFields()[cellI];
    GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];
    GaussField<tensor>& gradUG  = gradU.gaussFields()[cellI];
    GaussField<vector>& gradTG  = gradT.gaussFields()[cellI];

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
    // Convective term assembly:
    // rhoUG already stores rho*U, so no auxiliary flux field is needed.
    //
    // Diffusive term assembly:
    // no Navier-Stokes-Fourier mass diffusion term is used.
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
    // Build Fc_m at interior, minus, and plus Gauss points.
    //
    // Explicitly fill cell, minus, and plus Gauss values instead of using
    // dgExpr. This keeps exactly one materialized flux field.
    tmp<GaussField<tensor>> tMomentumFlux =
        GaussField<tensor>::New(cellI, &dgMesh);
    GaussField<tensor>& momentumFlux = tMomentumFlux.ref();

    for (label gpI = 0; gpI < momentumFlux.cellField().size(); ++gpI)
    {
        momentumFlux.cellField()[gpI] =
            rhoUG.cellField()[gpI]*UG.cellField()[gpI]
          + pG.cellField()[gpI]*symmTensor::I;
    }

    faceGaussField<tensor>& momentumFluxFace = momentumFlux.faceField();
    const faceGaussField<vector>& UFace = UG.faceField();
    const faceGaussField<vector>& rhoUFace = rhoUG.faceField();
    const faceGaussField<scalar>& pFace = pG.faceField();

    for (label faceI = 0; faceI < momentumFluxFace.nFaces(); ++faceI)
    {
        for (label gpI = 0; gpI < momentumFluxFace.nGaussPerFace(); ++gpI)
        {
            momentumFluxFace.minusValueOnFace(faceI, gpI) =
                rhoUFace.minusValueOnFace(faceI, gpI)
               *UFace.minusValueOnFace(faceI, gpI)
              + pFace.minusValueOnFace(faceI, gpI)*symmTensor::I;

            momentumFluxFace.plusValueOnFace(faceI, gpI) =
                rhoUFace.plusValueOnFace(faceI, gpI)
               *UFace.plusValueOnFace(faceI, gpI)
              + pFace.plusValueOnFace(faceI, gpI)*symmTensor::I;
        }
    }

    dgConvectiveTerm<vector, tensor> momentumConvTerm
    (
        "momentumConvTerm",
        dgFluxSolver::equationType::momentumTransport,
        basisField,
        rhoUG,
        tMomentumFlux,
        fluxSolverManager
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
    // Build Fc_E at interior, minus, and plus Gauss points.
    tmp<GaussField<vector>> tEnergyFlux =
        GaussField<vector>::New(cellI, &dgMesh);
    GaussField<vector>& energyFlux = tEnergyFlux.ref();

    for (label gpI = 0; gpI < energyFlux.cellField().size(); ++gpI)
    {
        energyFlux.cellField()[gpI] =
            (EG.cellField()[gpI] + pG.cellField()[gpI])
           *UG.cellField()[gpI];
    }

    faceGaussField<vector>& energyFluxFace = energyFlux.faceField();
    const faceGaussField<scalar>& EFace = EG.faceField();

    for (label faceI = 0; faceI < energyFluxFace.nFaces(); ++faceI)
    {
        for (label gpI = 0; gpI < energyFluxFace.nGaussPerFace(); ++gpI)
        {
            energyFluxFace.minusValueOnFace(faceI, gpI) =
                (
                    EFace.minusValueOnFace(faceI, gpI)
                  + pFace.minusValueOnFace(faceI, gpI)
                )
               *UFace.minusValueOnFace(faceI, gpI);

            energyFluxFace.plusValueOnFace(faceI, gpI) =
                (
                    EFace.plusValueOnFace(faceI, gpI)
                  + pFace.plusValueOnFace(faceI, gpI)
                )
               *UFace.plusValueOnFace(faceI, gpI);
        }
    }

    dgConvectiveTerm<scalar, vector> energyConvTerm
    (
        "energyConvTerm",
        dgFluxSolver::equationType::energyTransport,
        basisField,
        EG,
        tEnergyFlux,
        fluxSolverManager
    );

    //==================================================================//
    //             Shared Navier-Stokes-Fourier diffusive closure
    //==================================================================//
    //
    // The momentum and energy diffusive sub-blocks share tau.  The optimized
    // path computes both diffusive flux fields in one pass over Gauss points:
    //
    //     Fd_m = -tau,
    //     Fd_E = -(tau.U + kappa grad(T)).
    //
    // Transport properties from the thermo API:
    //     mu = dynamic viscosity
    //     kappa = thermal conductivity
    GaussField<scalar> muG(cellI, &dgMesh);
    GaussField<scalar> kappaG(cellI, &dgMesh);

    thermo->transport().calcMu(cellI, TG, muG);
    thermo->transport().calcKappa(cellI, TG, kappaG);

    // Momentum diffusive term assembly:
    // tMomentumDiffFlux stores Fd_m = -tau.
    //
    // Energy diffusive term assembly:
    // tEnergyDiffFlux stores Fd_E = -(tau.U + kappa grad(T)).
    tmp<GaussField<tensor>> tMomentumDiffFlux =
        GaussField<tensor>::New(cellI, &dgMesh);
    GaussField<tensor>& momentumDiffFlux = tMomentumDiffFlux.ref();

    tmp<GaussField<vector>> tEnergyDiffFlux =
        GaussField<vector>::New(cellI, &dgMesh);
    GaussField<vector>& energyDiffFlux = tEnergyDiffFlux.ref();

    for (label gpI = 0; gpI < momentumDiffFlux.cellField().size(); ++gpI)
    {
        const scalar muValue = muG.cellField()[gpI];
        const tensor& gradUValue = gradUG.cellField()[gpI];
        const scalar divU = tr(gradUValue);

        tensor tau = muValue*(gradUValue + Foam::T(gradUValue));
        tau.xx() -= scalar(2.0/3.0)*muValue*divU;
        tau.yy() -= scalar(2.0/3.0)*muValue*divU;
        tau.zz() -= scalar(2.0/3.0)*muValue*divU;

        momentumDiffFlux.cellField()[gpI] = -tau;
        energyDiffFlux.cellField()[gpI] =
           -(
                (tau & UG.cellField()[gpI])
              + kappaG.cellField()[gpI]*gradTG.cellField()[gpI]
            );
    }

    faceGaussField<tensor>& momentumDiffFace = momentumDiffFlux.faceField();
    faceGaussField<vector>& energyDiffFace = energyDiffFlux.faceField();

    const faceGaussField<scalar>& muFace = muG.faceField();
    const faceGaussField<scalar>& kappaFace = kappaG.faceField();
    const faceGaussField<tensor>& gradUFace = gradUG.faceField();
    const faceGaussField<vector>& gradTFace = gradTG.faceField();

    for (label faceI = 0; faceI < momentumDiffFace.nFaces(); ++faceI)
    {
        for (label gpI = 0; gpI < momentumDiffFace.nGaussPerFace(); ++gpI)
        {
            const scalar muMinus = muFace.minusValueOnFace(faceI, gpI);
            const tensor& gradUMinus =
                gradUFace.minusValueOnFace(faceI, gpI);
            const scalar divUMinus = tr(gradUMinus);

            tensor tauMinus = muMinus*(gradUMinus + Foam::T(gradUMinus));
            tauMinus.xx() -= scalar(2.0/3.0)*muMinus*divUMinus;
            tauMinus.yy() -= scalar(2.0/3.0)*muMinus*divUMinus;
            tauMinus.zz() -= scalar(2.0/3.0)*muMinus*divUMinus;

            const scalar muPlus = muFace.plusValueOnFace(faceI, gpI);
            const tensor& gradUPlus =
                gradUFace.plusValueOnFace(faceI, gpI);
            const scalar divUPlus = tr(gradUPlus);

            tensor tauPlus = muPlus*(gradUPlus + Foam::T(gradUPlus));
            tauPlus.xx() -= scalar(2.0/3.0)*muPlus*divUPlus;
            tauPlus.yy() -= scalar(2.0/3.0)*muPlus*divUPlus;
            tauPlus.zz() -= scalar(2.0/3.0)*muPlus*divUPlus;

            momentumDiffFace.minusValueOnFace(faceI, gpI) = -tauMinus;
            momentumDiffFace.plusValueOnFace(faceI, gpI) = -tauPlus;

            energyDiffFace.minusValueOnFace(faceI, gpI) =
               -(
                    (tauMinus & UFace.minusValueOnFace(faceI, gpI))
                  + kappaFace.minusValueOnFace(faceI, gpI)
                   *gradTFace.minusValueOnFace(faceI, gpI)
                );

            energyDiffFace.plusValueOnFace(faceI, gpI) =
               -(
                    (tauPlus & UFace.plusValueOnFace(faceI, gpI))
                  + kappaFace.plusValueOnFace(faceI, gpI)
                   *gradTFace.plusValueOnFace(faceI, gpI)
                );
        }
    }

    // Momentum diffusive DG term:
    // Integrate Fd_m with dgDiffusiveTerm's central face flux.
    dgDiffusiveTerm<vector, tensor> momentumDiffTerm
    (
        "momentumDiffTerm",
        dgFluxSolver::equationType::momentumTransport,
        basisField,
        rhoUG,
        tMomentumDiffFlux
    );

    // Energy diffusive DG term:
    // Integrate Fd_E with dgDiffusiveTerm's central face flux.
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
    // Numerical face fluxes:
    // - convective terms use dgFluxSolverManager-selected Riemann solvers,
    // - diffusive terms use dgDiffusiveTerm's fixed central flux.
    massConvTerm.updateFlux();
    massConvTerm.calcResidual();

    momentumConvTerm.updateFlux();
    momentumDiffTerm.updateFlux();
    momentumConvTerm.calcResidual();
    momentumDiffTerm.calcResidual();

    energyConvTerm.updateFlux();
    energyDiffTerm.updateFlux();
    energyConvTerm.calcResidual();
    energyDiffTerm.calcResidual();

    //==================================================================//
    //                   Store residuals for time update
    //==================================================================//
    //
    // Store each equation contribution on the matching conservative field.
    // Momentum and energy each receive both convective and diffusive residuals.
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
