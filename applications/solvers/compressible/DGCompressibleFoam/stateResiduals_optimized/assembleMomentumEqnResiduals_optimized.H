//======================================================================//
//                      Momentum conservation
//======================================================================//
//
//     d(rhoU)/dt + div(rhoU (x) U + pI + tauT) = rho*g.
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

tmp<GaussField<tensor>> tMomentumDiffFlux =
    GaussField<tensor>::New(cellI, &dgMesh);
GaussField<tensor>& momentumDiffFlux = tMomentumDiffFlux.ref();

tmp<GaussField<vector>> tMomentumSource =
    GaussField<vector>::New(cellI, &dgMesh);
GaussField<vector>& momentumSource = tMomentumSource.ref();

for (label gpI = 0; gpI < momentumDiffFlux.cellField().size(); ++gpI)
{
    const scalar muValue = muG.cellField()[gpI];
    const tensor& gradUValue = gradUG.cellField()[gpI];
    const scalar divU = tr(gradUValue);

    tensor tau = muValue*(gradUValue + Foam::T(gradUValue));
    tau.xx() -= scalar(2.0/3.0)*muValue*divU;
    tau.yy() -= scalar(2.0/3.0)*muValue*divU;
    tau.zz() -= scalar(2.0/3.0)*muValue*divU;

    const vector& mD = massDiffFlux.cellField()[gpI];
    const vector& UValue = UG.cellField()[gpI];
    tensor massMomentumDiff = mD*UValue + UValue*mD;
    const scalar mDDotU = mD & UValue;
    massMomentumDiff.xx() -= scalar(2.0/3.0)*mDDotU;
    massMomentumDiff.yy() -= scalar(2.0/3.0)*mDDotU;
    massMomentumDiff.zz() -= scalar(2.0/3.0)*mDDotU;

    momentumDiffFlux.cellField()[gpI] = -tau + massMomentumDiff;
    momentumSource.cellField()[gpI] = rhoG.cellField()[gpI]*g.value();
}

faceGaussField<tensor>& momentumDiffFace = momentumDiffFlux.faceField();
const faceGaussField<tensor>& gradUFace = gradUG.faceField();

for (label faceI = 0; faceI < momentumDiffFace.nFaces(); ++faceI)
{
    for (label gpI = 0; gpI < momentumDiffFace.nGaussPerFace(); ++gpI)
    {
        const scalar muMinus = muFace.minusValueOnFace(faceI, gpI);
        const tensor& gradUMinus = gradUFace.minusValueOnFace(faceI, gpI);
        const scalar divUMinus = tr(gradUMinus);

        tensor tauMinus = muMinus*(gradUMinus + Foam::T(gradUMinus));
        tauMinus.xx() -= scalar(2.0/3.0)*muMinus*divUMinus;
        tauMinus.yy() -= scalar(2.0/3.0)*muMinus*divUMinus;
        tauMinus.zz() -= scalar(2.0/3.0)*muMinus*divUMinus;

        const scalar muPlus = muFace.plusValueOnFace(faceI, gpI);
        const tensor& gradUPlus = gradUFace.plusValueOnFace(faceI, gpI);
        const scalar divUPlus = tr(gradUPlus);

        tensor tauPlus = muPlus*(gradUPlus + Foam::T(gradUPlus));
        tauPlus.xx() -= scalar(2.0/3.0)*muPlus*divUPlus;
        tauPlus.yy() -= scalar(2.0/3.0)*muPlus*divUPlus;
        tauPlus.zz() -= scalar(2.0/3.0)*muPlus*divUPlus;

        const vector& mDMinus =
            massDiffFace.minusValueOnFace(faceI, gpI);
        const vector& UMinus = UFace.minusValueOnFace(faceI, gpI);
        tensor massMomentumDiffMinus = mDMinus*UMinus + UMinus*mDMinus;
        const scalar mDMinusDotU = mDMinus & UMinus;
        massMomentumDiffMinus.xx() -= scalar(2.0/3.0)*mDMinusDotU;
        massMomentumDiffMinus.yy() -= scalar(2.0/3.0)*mDMinusDotU;
        massMomentumDiffMinus.zz() -= scalar(2.0/3.0)*mDMinusDotU;

        const vector& mDPlus = massDiffFace.plusValueOnFace(faceI, gpI);
        const vector& UPlus = UFace.plusValueOnFace(faceI, gpI);
        tensor massMomentumDiffPlus = mDPlus*UPlus + UPlus*mDPlus;
        const scalar mDPlusDotU = mDPlus & UPlus;
        massMomentumDiffPlus.xx() -= scalar(2.0/3.0)*mDPlusDotU;
        massMomentumDiffPlus.yy() -= scalar(2.0/3.0)*mDPlusDotU;
        massMomentumDiffPlus.zz() -= scalar(2.0/3.0)*mDPlusDotU;

        momentumDiffFace.minusValueOnFace(faceI, gpI) =
            -tauMinus + massMomentumDiffMinus;
        momentumDiffFace.plusValueOnFace(faceI, gpI) =
            -tauPlus + massMomentumDiffPlus;
    }
}

dgDiffusiveTerm<vector, tensor> momentumDiffTerm
(
    "momentumDiffTerm",
    dgFluxSolver::equationType::momentumTransport,
    basisField,
    rhoUG,
    tMomentumDiffFlux
);

dgSourceTerm<vector> momentumSourceTerm
(
    "momentumSourceTerm",
    basisField,
    momentumSource
);

momentumConvTerm.updateFlux();
momentumDiffTerm.updateFlux();
momentumConvTerm.calcResidual();
momentumDiffTerm.calcResidual();
momentumSourceTerm.calcResidual();

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
    rhoU.name(),
    cellI,
    momentumSourceTerm.R()
);
