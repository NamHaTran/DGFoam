/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgCompressibleBoundaryManager.H"
#include "dgThermoConservative.H"
#include "dgField.H"
#include "GaussField.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleBoundaryManager, 0);

dgCompressibleBoundaryManager::dgCompressibleBoundaryManager
(
    const dictionary& boundaryDict,
    const dgThermoConservative& thermo,
    const dgGeomMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            "dgCompressibleBoundaryManager",
            mesh.getFvMesh().time().constant(),
            mesh.getFvMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    thermo_(thermo),
    patchToCondition_(mesh.getFvMesh().boundary().size(), -1),
    internalP_(Zero),
    internalT_(Zero),
    internalU_(Zero)
{
    readInternalField(boundaryDict);
    constructBoundaryConditions(boundaryDict);
    resolveCouplings();
}


dgCompressibleBoundaryManager::dgCompressibleBoundaryManager
(
    const IOobject& io,
    const dgThermoConservative& thermo,
    const dgGeomMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            "dgCompressibleBoundaryManager",
            mesh.getFvMesh().time().constant(),
            mesh.getFvMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    thermo_(thermo),
    patchToCondition_(mesh.getFvMesh().boundary().size(), -1),
    internalP_(Zero),
    internalT_(Zero),
    internalU_(Zero)
{
    IOdictionary boundaryDict(io);
    readInternalField(boundaryDict);
    constructBoundaryConditions(boundaryDict);
    resolveCouplings();
}


dgCompressibleBoundaryManager::dgCompressibleBoundaryManager
(
    const dictionary& boundaryDict,
    const dgGeomMesh& mesh
)
:
    dgCompressibleBoundaryManager
    (
        boundaryDict,
        mesh.getFvMesh().lookupObject<dgThermoConservative>
        (
            "dgThermoConservative"
        ),
        mesh
    )
{}


dgCompressibleBoundaryManager::dgCompressibleBoundaryManager
(
    const IOobject& io,
    const dgGeomMesh& mesh
)
:
    dgCompressibleBoundaryManager
    (
        io,
        mesh.getFvMesh().lookupObject<dgThermoConservative>
        (
            "dgThermoConservative"
        ),
        mesh
    )
{}


void dgCompressibleBoundaryManager::readInternalField(const dictionary& dict)
{
    if (!dict.found("internalField"))
    {
        FatalIOErrorInFunction(dict)
            << "Missing 'internalField' entry in boundaryConditions dictionary"
            << exit(FatalIOError);
    }

    const dictionary& internalDict = dict.subDict("internalField");

    internalP_ = internalDict.get<scalar>("p");
    internalT_ = internalDict.get<scalar>("T");
    internalU_ = internalDict.get<vector>("U");
}


void dgCompressibleBoundaryManager::constructBoundaryConditions
(
    const dictionary& dict
)
{
    const dictionary& bfDict = dict.subDict("boundaryField");
    const fvBoundaryMesh& patches = mesh_.getFvMesh().boundary();
    const label nPatches = patches.size();

    for (label patchI = 0; patchI < nPatches; ++patchI)
    {
        const fvPatch& patch = patches[patchI];
        const word& patchName = patch.name();

        if (dgBoundaryHelper::isGeneralPatch(patch))
        {
            if (!bfDict.found(patchName))
            {
                FatalIOErrorInFunction(dict)
                    << "BoundaryField entry missing for patch: "
                    << patchName << exit(FatalIOError);
            }

            const dictionary& patchDict = bfDict.subDict(patchName);
            const label conditionI = bConditions_.size();

            bConditions_.append
            (
                dgCompressibleBoundaryField::New
                (
                    patch,
                    mesh_,
                    thermo_,
                    patchDict
                )
            );

            patchToCondition_[patchI] = conditionI;
            bConditions_.last()->checkPatchType();
        }
    }

    bConditions_.shrink();
}


void dgCompressibleBoundaryManager::resolveCouplings()
{
    forAll(bConditions_, conditionI)
    {
        bConditions_[conditionI]->resolveCoupling(*this);
    }
}


bool dgCompressibleBoundaryManager::writeData(Ostream& os) const
{
    os  << "nBoundaryConditions " << bConditions_.size() << token::END_STATEMENT
        << nl;

    return os.good();
}


label dgCompressibleBoundaryManager::conditionIndex
(
    const word& patchName
) const
{
    const fvBoundaryMesh& patches = mesh_.getFvMesh().boundary();
    const label patchID = patches.findPatchID(patchName);

    if (patchID < 0)
    {
        return -1;
    }

    return patchToCondition_[patchID];
}


bool dgCompressibleBoundaryManager::foundBoundaryField
(
    const word& patchName
) const
{
    return conditionIndex(patchName) >= 0;
}


const dgCompressibleBoundaryField& dgCompressibleBoundaryManager::boundaryField
(
    const word& patchName
) const
{
    const label conditionI = conditionIndex(patchName);

    if (conditionI < 0)
    {
        FatalErrorInFunction
            << "No compressible DG boundary condition found for patch "
            << patchName << exit(FatalError);
    }

    return bConditions_[conditionI]();
}


dgCompressibleBoundaryField& dgCompressibleBoundaryManager::boundaryFieldRef
(
    const word& patchName
)
{
    const label conditionI = conditionIndex(patchName);

    if (conditionI < 0)
    {
        FatalErrorInFunction
            << "No compressible DG boundary condition found for patch "
            << patchName << exit(FatalError);
    }

    return bConditions_[conditionI]();
}


void dgCompressibleBoundaryManager::initializeConservatives
(
    dgField<scalar>& rho,
    dgField<vector>& rhoU,
    dgField<scalar>& E
) const
{
    const scalar rho0 = thermo_.eos().calcRhoFromPT(internalP_, internalT_);
    const scalar e0   = thermo_.calcHeFromRhoT(rho0, internalT_);

    rho  = rho0;
    rhoU = rho0*internalU_;
    E    = rho0*(e0 + 0.5*magSqr(internalU_));
}


void dgCompressibleBoundaryManager::update
(
    GaussField<scalar>& rhoG,
    GaussField<vector>& rhoUG,
    GaussField<scalar>& EG,
    GaussField<vector>& UG
) const
{
    const label cellID = rhoG.cellID();

    faceGaussField<scalar>& rhoFace = rhoG.faceField();
    faceGaussField<vector>& rhoUFace = rhoUG.faceField();
    faceGaussField<scalar>& EFace = EG.faceField();
    faceGaussField<vector>& UFace = UG.faceField();

    for (label faceI = 0; faceI < rhoFace.nFaces(); ++faceI)
    {
        if (!rhoFace.isBoundary(faceI))
        {
            continue;
        }

        const label patchID = rhoFace.getOwnerPatchID(faceI);
        const label conditionI = patchToCondition_[patchID];

        if (conditionI < 0)
        {
            continue;
        }

        const vector n = rhoFace.faces()[faceI]->normal();

        for (label g = 0; g < rhoFace.nGaussPerFace(); ++g)
        {
            bConditions_[conditionI]->updateConservativeGhostState
            (
                cellID,
                faceI,
                g,
                n,
                rhoFace.minusValueOnFace(faceI, g),
                rhoUFace.minusValueOnFace(faceI, g),
                EFace.minusValueOnFace(faceI, g),
                rhoFace.plusValueOnFace(faceI, g),
                rhoUFace.plusValueOnFace(faceI, g),
                EFace.plusValueOnFace(faceI, g)
            );

            UFace.plusValueOnFace(faceI, g) =
                rhoUFace.plusValueOnFace(faceI, g)
               /rhoFace.plusValueOnFace(faceI, g);

            bConditions_[conditionI]->updatePrimitiveBCValue
            (
                cellID,
                faceI,
                g,
                n,
                rhoFace.minusValueOnFace(faceI, g),
                rhoUFace.minusValueOnFace(faceI, g),
                EFace.minusValueOnFace(faceI, g),
                rhoFace.midValueOnFace(faceI, g),
                rhoUFace.midValueOnFace(faceI, g),
                EFace.midValueOnFace(faceI, g)
            );
        }
    }
}


void dgCompressibleBoundaryManager::correctFlux
(
    GaussField<vector>& flux
) const
{
    const label cellID = flux.cellID();
    faceGaussField<vector>& fluxFace = flux.faceField();

    for (label faceI = 0; faceI < fluxFace.nFaces(); ++faceI)
    {
        if (!fluxFace.isBoundary(faceI))
        {
            continue;
        }

        const label patchID = fluxFace.getOwnerPatchID(faceI);
        const label conditionI = patchToCondition_[patchID];

        if (conditionI < 0)
        {
            continue;
        }

        const vector n = fluxFace.faces()[faceI]->normal();

        for (label g = 0; g < fluxFace.nGaussPerFace(); ++g)
        {
            bConditions_[conditionI]->correctFlux
            (
                cellID,
                faceI,
                g,
                n,
                fluxFace.minusValueOnFace(faceI, g)
            );

            fluxFace.plusValueOnFace(faceI, g) =
                fluxFace.minusValueOnFace(faceI, g);
        }
    }
}


void dgCompressibleBoundaryManager::correctGradient
(
    GaussField<vector>& grad
) const
{
    const label cellID = grad.cellID();
    faceGaussField<vector>& gradFace = grad.faceField();

    for (label faceI = 0; faceI < gradFace.nFaces(); ++faceI)
    {
        if (!gradFace.isBoundary(faceI))
        {
            continue;
        }

        const label patchID = gradFace.getOwnerPatchID(faceI);
        const label conditionI = patchToCondition_[patchID];

        if (conditionI < 0)
        {
            continue;
        }

        const vector n = gradFace.faces()[faceI]->normal();

        for (label g = 0; g < gradFace.nGaussPerFace(); ++g)
        {
            bConditions_[conditionI]->correctGradient
            (
                cellID,
                faceI,
                g,
                n,
                gradFace.minusValueOnFace(faceI, g)
            );

            gradFace.plusValueOnFace(faceI, g) =
                gradFace.minusValueOnFace(faceI, g);
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
