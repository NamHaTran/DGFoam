//======================================================================//
//                    Total-energy conservation
//======================================================================//
//
//     d(rho Ec)/dt + div(rho U Ec)
//       = div(qT) + Ek div(mD) - div(P U) - div(tauT.U).
//
// This explicit-loop version mirrors stateResiduals/assembleEnergyEqnResiduals.H.
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

for (label gpI = 0; gpI < energyDiffFlux.cellField().size(); ++gpI)
{
    energyDiffFlux.cellField()[gpI] =
        (momentumDiffFlux.cellField()[gpI] & UG.cellField()[gpI])
      - kappaG.cellField()[gpI]*gradTG.cellField()[gpI];
}

faceGaussField<vector>& energyDiffFace = energyDiffFlux.faceField();
const faceGaussField<scalar>& kappaFace = kappaG.faceField();

for (label faceI = 0; faceI < energyDiffFace.nFaces(); ++faceI)
{
    for (label gpI = 0; gpI < energyDiffFace.nGaussPerFace(); ++gpI)
    {
        energyDiffFace.minusValueOnFace(faceI, gpI) =
            (
                momentumDiffFace.minusValueOnFace(faceI, gpI)
              & UFace.minusValueOnFace(faceI, gpI)
            )
          - kappaFace.minusValueOnFace(faceI, gpI)
           *gradTFace.minusValueOnFace(faceI, gpI);

        energyDiffFace.plusValueOnFace(faceI, gpI) =
            (
                momentumDiffFace.plusValueOnFace(faceI, gpI)
              & UFace.plusValueOnFace(faceI, gpI)
            )
          - kappaFace.plusValueOnFace(faceI, gpI)
           *gradTFace.plusValueOnFace(faceI, gpI);
    }
}

tmp<GaussField<vector>> tSelfDiffusionEnergyFlux =
    GaussField<vector>::New(cellI, &dgMesh);
GaussField<vector>& selfDiffusionEnergyFlux =
    tSelfDiffusionEnergyFlux.ref();

if (selfDiffusion)
{
    for
    (
        label gpI = 0;
        gpI < selfDiffusionEnergyFlux.cellField().size();
        ++gpI
    )
    {
        const vector& UValue = UG.cellField()[gpI];
        const vector& mD = massDiffFlux.cellField()[gpI];
        const scalar internalEnergy =
            thermo->energyModel().calcHe(TG.cellField()[gpI]);

        selfDiffusionEnergyFlux.cellField()[gpI] =
            (selfDiffusionMomentumFlux.cellField()[gpI] & UValue)
          + internalEnergy*mD;
    }

    faceGaussField<vector>& selfDiffusionEnergyFace =
        selfDiffusionEnergyFlux.faceField();
    faceGaussField<tensor>& selfDiffusionMomentumFace =
        selfDiffusionMomentumFlux.faceField();

    for
    (
        label faceI = 0;
        faceI < selfDiffusionEnergyFace.nFaces();
        ++faceI
    )
    {
        for
        (
            label gpI = 0;
            gpI < selfDiffusionEnergyFace.nGaussPerFace();
            ++gpI
        )
        {
            selfDiffusionEnergyFace.minusValueOnFace(faceI, gpI) =
                (
                    selfDiffusionMomentumFace.minusValueOnFace(faceI, gpI)
                  & UFace.minusValueOnFace(faceI, gpI)
                )
              + (
                    thermo->energyModel().calcHe
                    (
                        TFace.minusValueOnFace(faceI, gpI)
                    )
                )*massDiffFace.minusValueOnFace(faceI, gpI);

            selfDiffusionEnergyFace.plusValueOnFace(faceI, gpI) =
                (
                    selfDiffusionMomentumFace.plusValueOnFace(faceI, gpI)
                  & UFace.plusValueOnFace(faceI, gpI)
                )
              + (
                    thermo->energyModel().calcHe
                    (
                        TFace.plusValueOnFace(faceI, gpI)
                    )
                )*massDiffFace.plusValueOnFace(faceI, gpI);
        }
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
        fluxSolverManager
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
        fluxSolverManager
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
    tmp<GaussField<vector>> tArtificialDiffusionEnergyFlux =
        GaussField<vector>::New(cellI, &dgMesh);
    GaussField<vector>& artificialDiffusionEnergyFlux =
        tArtificialDiffusionEnergyFlux.ref();

    for
    (
        label gpI = 0;
        gpI < artificialDiffusionEnergyFlux.cellField().size();
        ++gpI
    )
    {
        artificialDiffusionEnergyFlux.cellField()[gpI] =
           -muAVG.cellField()[gpI]*SEG.cellField()[gpI];
    }

    faceGaussField<vector>& artificialDiffusionEnergyFace =
        artificialDiffusionEnergyFlux.faceField();
    const faceGaussField<scalar>& muAVFace = muAVG.faceField();
    const faceGaussField<vector>& SEFace = SEG.faceField();

    for
    (
        label faceI = 0;
        faceI < artificialDiffusionEnergyFace.nFaces();
        ++faceI
    )
    {
        for
        (
            label gpI = 0;
            gpI < artificialDiffusionEnergyFace.nGaussPerFace();
            ++gpI
        )
        {
            artificialDiffusionEnergyFace.minusValueOnFace(faceI, gpI) =
               -muAVFace.minusValueOnFace(faceI, gpI)
               *SEFace.minusValueOnFace(faceI, gpI);

            artificialDiffusionEnergyFace.plusValueOnFace(faceI, gpI) =
               -muAVFace.plusValueOnFace(faceI, gpI)
               *SEFace.plusValueOnFace(faceI, gpI);
        }
    }

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
