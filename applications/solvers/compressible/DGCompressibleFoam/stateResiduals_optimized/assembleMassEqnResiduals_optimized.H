//======================================================================//
//                         Mass conservation
//======================================================================//
//
//     d(rho)/dt + div(Fc_rho + Fd_rho) = 0,
//     mD = -D0*mu*(grad(p)/p - 0.5*grad(T)/T).
dgConvectiveTerm<scalar, vector> massConvTerm
(
    "massConvTerm",
    dgFluxSolver::equationType::massTransport,
    dgFluxSolver::fluxType::convective,
    basisField,
    rhoG,
    rhoUG,
    fluxSolverManager
);

tmp<GaussField<vector>> tMassDiffFlux =
    GaussField<vector>::New(cellI, &dgMesh);
GaussField<vector>& massDiffFlux = tMassDiffFlux.ref();

faceGaussField<vector>& massDiffFace = massDiffFlux.faceField();
const faceGaussField<scalar>& muFace = muG.faceField();
const faceGaussField<vector>& gradPFace = gradPG.faceField();
const faceGaussField<vector>& gradTFace = gradTG.faceField();
const faceGaussField<scalar>& pFace = pG.faceField();
const faceGaussField<scalar>& TFace = TG.faceField();

massConvTerm.updateFlux();
massConvTerm.calcResidual();

if (selfDiffusion)
{
    for (label gpI = 0; gpI < massDiffFlux.cellField().size(); ++gpI)
    {
        massDiffFlux.cellField()[gpI] =
           -selfDiffCoeffScale*muG.cellField()[gpI]
           *(
                gradPG.cellField()[gpI]/pG.cellField()[gpI]
              - scalar(0.5)*gradTG.cellField()[gpI]/TG.cellField()[gpI]
            );
    }

    for (label faceI = 0; faceI < massDiffFace.nFaces(); ++faceI)
    {
        for (label gpI = 0; gpI < massDiffFace.nGaussPerFace(); ++gpI)
        {
            massDiffFace.minusValueOnFace(faceI, gpI) =
               -selfDiffCoeffScale*muFace.minusValueOnFace(faceI, gpI)
               *(
                    gradPFace.minusValueOnFace(faceI, gpI)
                   /pFace.minusValueOnFace(faceI, gpI)
                  - scalar(0.5)
                   *gradTFace.minusValueOnFace(faceI, gpI)
                   /TFace.minusValueOnFace(faceI, gpI)
                );

            massDiffFace.plusValueOnFace(faceI, gpI) =
               -selfDiffCoeffScale*muFace.plusValueOnFace(faceI, gpI)
               *(
                    gradPFace.plusValueOnFace(faceI, gpI)
                   /pFace.plusValueOnFace(faceI, gpI)
                  - scalar(0.5)
                   *gradTFace.plusValueOnFace(faceI, gpI)
                   /TFace.plusValueOnFace(faceI, gpI)
                );
        }
    }

    compressibleBC->correctFlux(massDiffFlux);

    dgDiffusiveTerm<scalar, vector> selfDiffusionMassTerm
    (
        "selfDiffusionMassTerm",
        dgFluxSolver::equationType::massTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        rhoG,
        massDiffFlux,
        fluxSolverManager
    );

    selfDiffusionMassTerm.updateFlux();
    selfDiffusionMassTerm.calcResidual();

    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        massConvTerm.R()
    );

    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        selfDiffusionMassTerm.R()
    );
}
else
{
    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        massConvTerm.R()
    );
}

if (artificialDiffusion)
{
    tmp<GaussField<vector>> tArtificialDiffusionMassFlux =
        GaussField<vector>::New(cellI, &dgMesh);
    GaussField<vector>& artificialDiffusionMassFlux =
        tArtificialDiffusionMassFlux.ref();

    for
    (
        label gpI = 0;
        gpI < artificialDiffusionMassFlux.cellField().size();
        ++gpI
    )
    {
        artificialDiffusionMassFlux.cellField()[gpI] =
           -muAVG.cellField()[gpI]*SRhoG.cellField()[gpI];
    }

    faceGaussField<vector>& artificialDiffusionMassFace =
        artificialDiffusionMassFlux.faceField();
    const faceGaussField<scalar>& muAVFace = muAVG.faceField();
    const faceGaussField<vector>& SRhoFace = SRhoG.faceField();

    for
    (
        label faceI = 0;
        faceI < artificialDiffusionMassFace.nFaces();
        ++faceI
    )
    {
        for
        (
            label gpI = 0;
            gpI < artificialDiffusionMassFace.nGaussPerFace();
            ++gpI
        )
        {
            artificialDiffusionMassFace.minusValueOnFace(faceI, gpI) =
               -muAVFace.minusValueOnFace(faceI, gpI)
               *SRhoFace.minusValueOnFace(faceI, gpI);

            artificialDiffusionMassFace.plusValueOnFace(faceI, gpI) =
               -muAVFace.plusValueOnFace(faceI, gpI)
               *SRhoFace.plusValueOnFace(faceI, gpI);
        }
    }

    dgDiffusiveTerm<scalar, vector> artificialDiffusionMassTerm
    (
        "artificialDiffusionMassTerm",
        dgFluxSolver::equationType::massTransport,
        dgFluxSolver::fluxType::diffusive,
        basisField,
        rhoG,
        tArtificialDiffusionMassFlux,
        fluxSolverManager
    );

    artificialDiffusionMassTerm.updateFlux();
    artificialDiffusionMassTerm.calcResidual();

    timeDiscretization.addResidual
    (
        rho.name(),
        cellI,
        artificialDiffusionMassTerm.R()
    );
}
