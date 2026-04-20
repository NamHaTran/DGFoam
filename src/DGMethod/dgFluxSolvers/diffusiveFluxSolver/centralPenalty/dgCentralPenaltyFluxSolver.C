/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

\*---------------------------------------------------------------------------*/

#include "dgCentralPenaltyFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{

defineTypeNameAndDebug(dgCentralPenaltyFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgCentralPenaltyFluxSolver, dictionary);


word dgCentralPenaltyFluxSolver::readMuName(const dictionary& dict)
{
    const dictionary* readDict = &dict;

    if (dict.found("fluxSolversCoeffs"))
    {
        const dictionary& coeffsDict = dict.subDict("fluxSolversCoeffs");
        const word coeffKey("centralPenaltyCoeffs");

        if (coeffsDict.found(coeffKey))
        {
            readDict = &coeffsDict.subDict(coeffKey);
        }
    }

    return readDict->lookupOrDefault<word>("mu", "mu");
}


dgCentralPenaltyFluxSolver::dgCentralPenaltyFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgFluxSolver(name, dict, mesh),
    mu_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>
        (
            readMuName(dict)
        )
    ),
    CIP_(1.0),
    muName_(readMuName(dict))
{
    fType_ = dgFluxSolver::fluxType::diffusive;
    read(dict);
}


void dgCentralPenaltyFluxSolver::read(const dictionary& dict)
{
    const dictionary* readDict = &dict;

    if (dict.found("fluxSolversCoeffs"))
    {
        const dictionary& coeffsDict = dict.subDict("fluxSolversCoeffs");
        const word coeffKey("centralPenaltyCoeffs");

        if (coeffsDict.found(coeffKey))
        {
            readDict = &coeffsDict.subDict(coeffKey);
        }
    }

    CIP_ = readDict->lookupOrDefault<scalar>("CIP", CIP_);
    const word requestedMuName = readDict->lookupOrDefault<word>("mu", muName_);

    if (requestedMuName != muName_)
    {
        FatalIOErrorInFunction(*readDict)
            << "Cannot change viscosity field name for centralPenalty after "
            << "construction. Constructed with mu = " << muName_
            << ", requested mu = " << requestedMuName << "."
            << nl << exit(FatalIOError);
    }
}


scalar dgCentralPenaltyFluxSolver::penaltyCoeff
(
    const label cellID,
    const label fI,
    const faceGaussField<scalar>& muFace
) const
{
    const dgGeomMesh* dgMesh = muFace.dgMesh();
    const label faceID = muFace.globalFaceID(fI);
    const fvMesh& fvMesh = dgMesh->getFvMesh();

    const label owner = fvMesh.faceOwner()[faceID];
    const scalar VO = fvMesh.V()[owner];
    scalar VN = VO;

    if (faceID < fvMesh.nInternalFaces())
    {
        const label neighbour = fvMesh.faceNeighbour()[faceID];
        VN = fvMesh.V()[neighbour];
    }

    const scalar faceArea = max(muFace.faceAreas()[fI], VSMALL);
    const scalar h = max(min(VO, VN)/faceArea, VSMALL);
    const scalar p = scalar(dgMesh->pOrder());

    scalar mu = scalar(0);
    for (label gI = 0; gI < muFace.nGaussPerFace(); ++gI)
    {
        mu += scalar(0.5)
            *(
                muFace.minusValueOnFace(fI, gI)
              + muFace.plusValueOnFace(fI, gI)
             );
    }
    mu /= max(muFace.nGaussPerFace(), label(1));

    return CIP_*mu*sqr(p + scalar(1))/h;
}


void dgCentralPenaltyFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<vector>& F,
    const faceGaussField<scalar>& U
)
{
    const faceGaussField<scalar>& muFace =
        mu_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];
        const scalar penalty = penaltyCoeff(cellID, fI, muFace);

        for (label gI = 0; gI < nGauss; ++gI)
        {
            if (F.isBoundary(fI) && !F.isProcessorPatch(fI))
            {
                F.plusValueOnFace(fI, gI) = F.minusValueOnFace(fI, gI);
            }

            const vector FL = F.minusValueOnFace(fI, gI);
            const vector FR = F.plusValueOnFace(fI, gI);
            const scalar UL = U.minusValueOnFace(fI, gI);
            const scalar UR = U.plusValueOnFace(fI, gI);

            const scalar fn =
                scalar(0.5)*((FL & n) + (FR & n))
              - penalty*(UR - UL);

            F.fluxOnFace(fI, gI) = fn*n;
        }
    }
}


void dgCentralPenaltyFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>& U
)
{
    const faceGaussField<scalar>& muFace =
        mu_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];
        const scalar penalty = penaltyCoeff(cellID, fI, muFace);

        for (label gI = 0; gI < nGauss; ++gI)
        {
            if (F.isBoundary(fI) && !F.isProcessorPatch(fI))
            {
                F.plusValueOnFace(fI, gI) = F.minusValueOnFace(fI, gI);
            }

            const tensor FL = F.minusValueOnFace(fI, gI);
            const tensor FR = F.plusValueOnFace(fI, gI);
            const vector UL = U.minusValueOnFace(fI, gI);
            const vector UR = U.plusValueOnFace(fI, gI);

            const vector fn =
                scalar(0.5)*((FL & n) + (FR & n))
              - penalty*(UR - UL);

            F.fluxOnFace(fI, gI) = fn;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
