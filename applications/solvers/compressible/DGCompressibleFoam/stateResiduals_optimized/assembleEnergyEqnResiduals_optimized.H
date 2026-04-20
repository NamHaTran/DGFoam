//======================================================================//
//                    Total-energy conservation
//======================================================================//
//
//     d(E)/dt + div((E + p)U + tauT.U - (E/rho)mD - kappa gradT) = S_E.
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
    dgFluxSolver::fluxType::convective,
    basisField,
    EG,
    tEnergyFlux,
    fluxSolverManager
);

tmp<GaussField<vector>> tEnergyDiffFlux =
    GaussField<vector>::New(cellI, &dgMesh);
GaussField<vector>& energyDiffFlux = tEnergyDiffFlux.ref();

tmp<GaussField<scalar>> tEnergySource =
    GaussField<scalar>::New(cellI, &dgMesh);
GaussField<scalar>& energySource = tEnergySource.ref();

for (label gpI = 0; gpI < energyDiffFlux.cellField().size(); ++gpI)
{
    const scalar rhoValue = rhoG.cellField()[gpI];
    const scalar EValue = EG.cellField()[gpI];
    const vector& UValue = UG.cellField()[gpI];
    const vector& mD = massDiffFlux.cellField()[gpI];
    const tensor& tauT = momentumDiffFlux.cellField()[gpI];

    energyDiffFlux.cellField()[gpI] =
        (tauT & UValue)
      - (EValue/rhoValue)*mD
      - kappaG.cellField()[gpI]*gradTG.cellField()[gpI];

    const vector gradSpecificEnergy =
        (
            rhoValue*SEG.cellField()[gpI]
          - EValue*SRhoG.cellField()[gpI]
        )/(rhoValue*rhoValue);

    energySource.cellField()[gpI] =
        pG.cellField()[gpI]*tr(gradUG.cellField()[gpI])
      + (tauT && gradUG.cellField()[gpI])
      - (mD & gradSpecificEnergy)
      + rhoValue*(g.value() & UValue);
}

faceGaussField<vector>& energyDiffFace = energyDiffFlux.faceField();
const faceGaussField<scalar>& rhoFace = rhoG.faceField();
const faceGaussField<scalar>& kappaFace = kappaG.faceField();

for (label faceI = 0; faceI < energyDiffFace.nFaces(); ++faceI)
{
    for (label gpI = 0; gpI < energyDiffFace.nGaussPerFace(); ++gpI)
    {
        const tensor& tauTMinus =
            momentumDiffFace.minusValueOnFace(faceI, gpI);
        const vector& UMinus = UFace.minusValueOnFace(faceI, gpI);
        const vector& mDMinus = massDiffFace.minusValueOnFace(faceI, gpI);

        energyDiffFace.minusValueOnFace(faceI, gpI) =
            (tauTMinus & UMinus)
          - (
                EFace.minusValueOnFace(faceI, gpI)
               /rhoFace.minusValueOnFace(faceI, gpI)
            )*mDMinus
          - kappaFace.minusValueOnFace(faceI, gpI)
           *gradTFace.minusValueOnFace(faceI, gpI);

        const tensor& tauTPlus =
            momentumDiffFace.plusValueOnFace(faceI, gpI);
        const vector& UPlus = UFace.plusValueOnFace(faceI, gpI);
        const vector& mDPlus = massDiffFace.plusValueOnFace(faceI, gpI);

        energyDiffFace.plusValueOnFace(faceI, gpI) =
            (tauTPlus & UPlus)
          - (
                EFace.plusValueOnFace(faceI, gpI)
               /rhoFace.plusValueOnFace(faceI, gpI)
            )*mDPlus
          - kappaFace.plusValueOnFace(faceI, gpI)
           *gradTFace.plusValueOnFace(faceI, gpI);
    }
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

dgSourceTerm<scalar> energySourceTerm
(
    "energySourceTerm",
    basisField,
    energySource
);

energyConvTerm.updateFlux();
energyDiffTerm.updateFlux();
energyConvTerm.calcResidual();
energyDiffTerm.calcResidual();
energySourceTerm.calcResidual();

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
    energySourceTerm.R()
);
