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

dgCompressibleBoundaryManager::dgCompressibleBoundaryManager
(
    const dictionary& boundaryDict,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    thermo_
    (
        mesh.getFvMesh().lookupObject<dgThermoConservative>("dgThermoConservative")
    ),
    patchToCondition_(mesh.getFvMesh().boundary().size(), -1),
    internalP_(Zero),
    internalT_(Zero),
    internalU_(Zero)
{
    readInternalField(boundaryDict);
    constructBoundaryConditions(boundaryDict);
}


dgCompressibleBoundaryManager::dgCompressibleBoundaryManager
(
    const IOobject& io,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    thermo_
    (
        mesh.getFvMesh().lookupObject<dgThermoConservative>("dgThermoConservative")
    ),
    patchToCondition_(mesh.getFvMesh().boundary().size(), -1),
    internalP_(Zero),
    internalT_(Zero),
    internalU_(Zero)
{
    IOdictionary boundaryDict(io);
    readInternalField(boundaryDict);
    constructBoundaryConditions(boundaryDict);
}


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


void dgCompressibleBoundaryManager::updateValue
(
    GaussField<scalar>& rhoG,
    GaussField<vector>& rhoUG,
    GaussField<scalar>& EG
) const
{
    const label cellID = rhoG.cellID();

    faceGaussField<scalar>& rhoFace = rhoG.faceField();
    faceGaussField<vector>& rhoUFace = rhoUG.faceField();
    faceGaussField<scalar>& EFace = EG.faceField();

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
            bConditions_[conditionI]->updateGhostState
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
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
