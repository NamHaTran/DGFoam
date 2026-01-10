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

#ifndef Foam_dgGeneralBoundaryManagerI_H
#define Foam_dgGeneralBoundaryManagerI_H

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::dgGeneralBoundaryManager<Type>::dgGeneralBoundaryManager
(
    const dictionary& fieldDict,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh)
{
    // --- Read internalField entry
    const entry* internalEntry =
        fieldDict.findEntry("internalField", keyType::LITERAL);

    if (!internalEntry)
    {
        FatalIOErrorInFunction(fieldDict)
            << "Missing 'internalField' entry in dictionary"
            << exit(FatalIOError);
    }

    ITstream& is = const_cast<entry*>(internalEntry)->stream();

    word fieldType;
    is >> fieldType;

    if (fieldType != "uniform")
    {
        FatalIOErrorInFunction(fieldDict)
            << "Only 'uniform' internalField is supported. Got: "
            << fieldType << exit(FatalIOError);
    }

    is >> internalValue_;

    // --- Read boundaryField
    const dictionary& bfDict = fieldDict.subDict("boundaryField");

    const fvBoundaryMesh& patches = mesh_.getFvMesh().boundary();
    const label nPatches = patches.size();

    for (label patchI = 0; patchI < nPatches; ++patchI)
    {
        const fvPatch& patch = patches[patchI];
        const word& patchName = patch.name();

        if (!bfDict.found(patchName))
        {
            FatalIOErrorInFunction(fieldDict)
                << "BoundaryField entry missing for patch: "
                << patchName << exit(FatalIOError);
        }

        if (dgBoundaryHelper::isGeneralPatch(patch))
        {
            const dictionary& patchDict = bfDict.subDict(patchName);

            bConditions_.append
            (
                dgGeneralBoundaryField<Type>::New
                (
                    patch,
                    mesh_,
                    patchDict
                )
            );

            // Check validity of type
            bConditions_.last()->checkPatchType();
        }
    }

    bConditions_.shrink();
}


template<class Type>
Foam::dgGeneralBoundaryManager<Type>::dgGeneralBoundaryManager
(
    const IOobject& io,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh)
{
    IOdictionary fieldDict(io);

    // --- Read internalField entry
    const entry* internalEntry =
        fieldDict.findEntry("internalField", keyType::LITERAL);

    if (!internalEntry)
    {
        FatalIOErrorInFunction(fieldDict)
            << "Missing 'internalField' entry in dictionary"
            << exit(FatalIOError);
    }

    ITstream& is = const_cast<entry*>(internalEntry)->stream();

    word fieldType;
    is >> fieldType;

    if (fieldType != "uniform")
    {
        FatalIOErrorInFunction(fieldDict)
            << "Only 'uniform' internalField is supported. Got: "
            << fieldType << exit(FatalIOError);
    }

    is >> internalValue_;

    // --- Read boundaryField
    const dictionary& bfDict = fieldDict.subDict("boundaryField");

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
                FatalIOErrorInFunction(fieldDict)
                    << "BoundaryField entry missing for patch: "
                    << patchName << exit(FatalIOError);
            }

            const dictionary& patchDict = bfDict.subDict(patchName);

            bConditions_.append
            (
                dgGeneralBoundaryField<Type>::New
                (
                    patch,
                    mesh_,
                    patchDict
                )
            );

            // Check validity of type
            bConditions_.last()->checkPatchType();
        }
    }

    bConditions_.shrink();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void dgGeneralBoundaryManager<Type>::updateValue
(
    GaussField<Type>& GField
)
{
    //- Cell ID of this GaussField
    const label cellID = GField.cellID();

    //- Reference to faceGaussField of this cell
    faceGaussField<Type>& faceField = GField.faceField();

    //- Loop over local faces of this cell
    for (label faceI = 0; faceI < faceField.nFaces(); ++faceI)
    {
        //- Check if boundary face
        if (faceField.isBoundary(faceI))
        {
            //- Get owner patch ID from faceGaussField
            const label patchID = faceField.getOwnerPatchID(faceI);

            //- Skip non-general patches
            if (!this->isGeneralPatch(patchID))
            {
                continue;
            }

            //- Face unit normal vector
            const vector n = faceField.faces()[faceI]->normal();

            //- Loop over Gauss points on this face
            for (label g = 0; g < faceField.nGaussPerFace(); ++g)
            {
                const Type& mVal = faceField.minusValueOnFace(faceI, g);
                Type&       pVal = faceField.plusValueOnFace(faceI, g);

                //- Call BC implementation
                bConditions_[patchID]->updateGhostStateValue
                (
                    cellID,
                    faceI,      // local face index
                    g,          // Gauss point index
                    n,          // face normal
                    mVal,       // interior state U-
                    pVal        // ghost state U+ (output)
                );
            }
        }
    }
}

template<class Type>
void dgGeneralBoundaryManager<Type>::updateGrad
(
    GaussField<Type>& GField
)
{
    //- Cell ID of this GaussField
    const label cellID = GField.cellID();

    //- Reference to faceGaussField of this cell
    faceGaussField<Type>& faceField = GField.faceField();

    //- Loop over local faces
    for (label faceI = 0; faceI < faceField.nFaces(); ++faceI)
    {
        if (faceField.isBoundary(faceI))
        {
            const label patchID = faceField.getOwnerPatchID(faceI);

            //- Skip non-general patches
            if (!this->isGeneralPatch(patchID))
            {
                continue;
            }

            const vector n = faceField.faces()[faceI].normal();

            //- Loop over Gauss points
            for (label g = 0; g < faceField.nGaussPerFace(); ++g)
            {
                const Type& mVal = faceField.minusValueOnFace(faceI, g);
                Type&       pVal = faceField.plusValueOnFace(faceI, g);

                //- Call gradient BC implementation
                bConditions_[patchID].updateGhostStateGrad
                (
                    cellID,
                    faceI,
                    g,
                    n,
                    mVal,   // input S-
                    pVal    // output S+
                );
            }
        }
    }
}

#endif
// ************************************************************************* //
