//======================================================================//
//                         Mass conservation
//======================================================================//
//
//     d(rho)/dt + div(rho U + mD) = 0,
//     mD = -mu*(grad(p)/p - 0.5*grad(T)/T).
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

for (label gpI = 0; gpI < massDiffFlux.cellField().size(); ++gpI)
{
    massDiffFlux.cellField()[gpI] =
       -muG.cellField()[gpI]
       *(
            gradPG.cellField()[gpI]/pG.cellField()[gpI]
          - scalar(0.5)*gradTG.cellField()[gpI]/TG.cellField()[gpI]
        );
}

faceGaussField<vector>& massDiffFace = massDiffFlux.faceField();
const faceGaussField<scalar>& muFace = muG.faceField();
const faceGaussField<vector>& gradPFace = gradPG.faceField();
const faceGaussField<vector>& gradTFace = gradTG.faceField();
const faceGaussField<scalar>& pFace = pG.faceField();
const faceGaussField<scalar>& TFace = TG.faceField();

for (label faceI = 0; faceI < massDiffFace.nFaces(); ++faceI)
{
    for (label gpI = 0; gpI < massDiffFace.nGaussPerFace(); ++gpI)
    {
        massDiffFace.minusValueOnFace(faceI, gpI) =
           -muFace.minusValueOnFace(faceI, gpI)
           *(
                gradPFace.minusValueOnFace(faceI, gpI)
               /pFace.minusValueOnFace(faceI, gpI)
              - scalar(0.5)
               *gradTFace.minusValueOnFace(faceI, gpI)
               /TFace.minusValueOnFace(faceI, gpI)
            );

        massDiffFace.plusValueOnFace(faceI, gpI) =
           -muFace.plusValueOnFace(faceI, gpI)
           *(
                gradPFace.plusValueOnFace(faceI, gpI)
               /pFace.plusValueOnFace(faceI, gpI)
              - scalar(0.5)
               *gradTFace.plusValueOnFace(faceI, gpI)
               /TFace.plusValueOnFace(faceI, gpI)
            );
    }
}

dgDiffusiveTerm<scalar, vector> massDiffTerm
(
    "massDiffTerm",
    dgFluxSolver::equationType::massTransport,
    dgFluxSolver::fluxType::diffusive,
    basisField,
    rhoG,
    tMassDiffFlux,
    fluxSolverManager
);

massConvTerm.updateFlux();
massDiffTerm.updateFlux();
massConvTerm.calcResidual();
massDiffTerm.calcResidual();

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
    massDiffTerm.R()
);
