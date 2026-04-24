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


const dictionary& dgCentralPenaltyFluxSolver::coeffDict
(
    const word& name,
    const dictionary& dict
)
{
    if (dict.found("fluxSolversCoeffs"))
    {
        const dictionary& coeffsDict = dict.subDict("fluxSolversCoeffs");
        const word termCoeffKey(name + "Coeffs");

        if (coeffsDict.found(termCoeffKey))
        {
            return coeffsDict.subDict(termCoeffKey);
        }

        if (coeffsDict.found("centralPenaltyCoeffs"))
        {
            return coeffsDict.subDict("centralPenaltyCoeffs");
        }
    }

    return dict;
}


word dgCentralPenaltyFluxSolver::readDiffCoefName
(
    const word& name,
    const dictionary& dict
)
{
    const dictionary& readDict = coeffDict(name, dict);

    if (!readDict.found("diffCoef"))
    {
        FatalIOErrorInFunction(readDict)
            << "Missing required entry 'diffCoef' for centralPenalty flux "
            << "solver '" << name << "'."
            << nl << exit(FatalIOError);
    }

    return readDict.get<word>("diffCoef");
}


dgCentralPenaltyFluxSolver::dgCentralPenaltyFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgFluxSolver(name, dict, mesh),
    diffCoef_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>
        (
            readDiffCoefName(name, dict)
        )
    ),
    CIP_(1.0),
    diffCoefName_(readDiffCoefName(name, dict))
{
    fType_ = dgFluxSolver::fluxType::diffusive;
    read(dict);
}


void dgCentralPenaltyFluxSolver::read(const dictionary& dict)
{
    const dictionary& readDict = coeffDict(name_, dict);

    CIP_ = readDict.lookupOrDefault<scalar>("CIP", CIP_);
    const word requestedDiffCoefName =
        readDict.lookupOrDefault<word>("diffCoef", diffCoefName_);

    if (requestedDiffCoefName != diffCoefName_)
    {
        FatalIOErrorInFunction(readDict)
            << "Cannot change diffusion-coefficient field name for "
            << "centralPenalty after construction. Constructed with "
            << "diffCoef = " << diffCoefName_
            << ", requested diffCoef = " << requestedDiffCoefName << "."
            << nl << exit(FatalIOError);
    }
}


scalar dgCentralPenaltyFluxSolver::penaltyCoeff
(
    const label cellID,
    const label fI,
    const faceGaussField<scalar>& diffCoefFace
) const
{
    const dgGeomMesh* dgMesh = diffCoefFace.dgMesh();
    const label faceID = diffCoefFace.globalFaceID(fI);
    const fvMesh& fvMesh = dgMesh->getFvMesh();

    const label owner = fvMesh.faceOwner()[faceID];
    const scalar VO = fvMesh.V()[owner];
    scalar VN = VO;

    if (faceID < fvMesh.nInternalFaces())
    {
        const label neighbour = fvMesh.faceNeighbour()[faceID];
        VN = fvMesh.V()[neighbour];
    }

    const scalar faceArea = max(diffCoefFace.faceAreas()[fI], VSMALL);
    const scalar h = max(min(VO, VN)/faceArea, VSMALL);
    const scalar p = scalar(dgMesh->pOrder());

    scalar diffCoef = scalar(0);
    for (label gI = 0; gI < diffCoefFace.nGaussPerFace(); ++gI)
    {
        diffCoef += scalar(0.5)
            *(
                diffCoefFace.minusValueOnFace(fI, gI)
              + diffCoefFace.plusValueOnFace(fI, gI)
             );
    }
    diffCoef /= max(diffCoefFace.nGaussPerFace(), label(1));

    return CIP_*diffCoef*sqr(p + scalar(1))/h;
}


void dgCentralPenaltyFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<vector>& F,
    const faceGaussField<scalar>& U
)
{
    const faceGaussField<scalar>& diffCoefFace =
        diffCoef_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];
        const scalar penalty = penaltyCoeff(cellID, fI, diffCoefFace);

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
    const faceGaussField<scalar>& diffCoefFace =
        diffCoef_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];
        const scalar penalty = penaltyCoeff(cellID, fI, diffCoefFace);

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
