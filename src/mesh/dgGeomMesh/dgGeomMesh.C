/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgGeomMesh.H"
#include "dgGeomFace.H"
#include "dgRefFace.H"
#include "dgGeomCell.H"

#include "faceList.H"
#include "cellList.H"

#include "polyBoundaryMesh.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from fvMesh and polynomial order
Foam::dgGeomMesh::dgGeomMesh
(
    const fvMesh& mesh,
    const label pOrder
)
:
    mesh_(mesh),
    pOrder_(pOrder),
    refFace_(std::make_shared<dgRefFace>(pOrder)),
    refCellTet_(std::make_shared<dgRefCell>(pOrder,dgCellType::TET)),
    refCellHex_(std::make_shared<dgRefCell>(pOrder,dgCellType::HEX)),
    refCellPrism_(std::make_shared<dgRefCell>(pOrder,dgCellType::PRISM)),
    refCellPyramid_(std::make_shared<dgRefCell>(pOrder,dgCellType::PYRAMID))
{
    // Create all geometric faces

    const label nF = mesh_.nFaces();
    faces_.setSize(nF);

    // Create dgGeomFace objects for each face
    for (label faceI = 0; faceI < nF; ++faceI)
    {
        // Allocate new dgGeomFace using shared refFace_
        faces_[faceI] = new dgGeomFace(faceI, mesh_, refFace_);
    }

    // Create dgGeomCell objects
    cells_.setSize(mesh_.nCells());
    for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
    {
        // Create dgGeomCell for each cell
        const labelList& cellPoints = mesh.cellPoints()[cellI];
        const label nPoints = cellPoints.size();

        // Guess the type of cell based on the number of points
        std::shared_ptr<dgRefCell> refCell = nullptr;

        if (nPoints == 4)
        {
            refCell = refCellTet_;
        }
        else if (nPoints == 8)
        {
            refCell = refCellHex_;
        }
        else if (nPoints == 6)
        {
            refCell = refCellPrism_;
        }
        else if (nPoints == 5)
        {
            refCell = refCellPyramid_;
        }
        else
        {
            Info << "Skip unsupported cell with " << nPoints << " points." << endl;
            continue;
        }

        // Create geomCell
        cells_[cellI] = new dgGeomCell(cellI, mesh_, refCell);
        // cells_[cellI]->printDebugInfo();
        cells_[cellI]->updateFaceInfo(faces_);
    }

    // Update face connectivity and calculate basis functions
    for (label faceI = 0; faceI < nF; ++faceI)
    {
        faces_[faceI]->findGaussConnectivity();

        // Note, computing basis functions and derivatives must be done after
        // the connectivity is established, as it depends on the Gauss points.
        faces_[faceI]->computeBasisAndDerivatives();
    }

    // Get boundary faces
    getBoundaryFaces();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::dgGeomMesh::~dgGeomMesh()
{
    // Delete all dynamically allocated face objects
    forAll(faces_, i)
    {
        delete faces_[i];
    }

    // Delete all dynamically allocated cell objects
    forAll(cells_, i)
    {
        delete cells_[i];
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dgGeomMesh::getBoundaryFaces()
{
    // Clear the list in case it was previously filled
    boundaryFaces_.clear();

    // Get the list of all patches in the mesh
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Temporary dynamic list to hold global face IDs on the boundary
    DynamicList<label> bFaces;

    // Loop through each patch
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Skip patches of type 'empty' (used in 2D) and 'processor' (used in parallel)
        if (isA<emptyPolyPatch>(pp) || isA<processorPolyPatch>(pp))
        {
            continue;
        }

        // Get the starting global face ID of this patch
        const label start = pp.start();

        // Get the number of faces in this patch
        const label size = pp.size();

        // Add each face ID in this patch to the boundary list
        for (label i = 0; i < size; ++i)
        {
            const label faceID = start + i;
            bFaces.append(faceID);
        }
    }

    // Transfer the result into the class member boundaryFaces_
    boundaryFaces_.transfer(bFaces);
}

const Foam::labelList& Foam::dgGeomMesh::boundaryFaces() const
{
    return boundaryFaces_;
}

Foam::label Foam::dgGeomMesh::getPatchID(const label faceID) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Skip 'empty' and 'processor' patch types
        if (isA<emptyPolyPatch>(pp) || isA<processorPolyPatch>(pp))
        {
            continue;
        }

        const label start = pp.start();
        const label size  = pp.size();

        // Check if faceID falls within this patch
        if (faceID >= start && faceID < start + size)
        {
            return patchI;
        }
    }

    // If not found, throw an error
    FatalErrorInFunction
        << "Face ID " << faceID << " does not belong to any valid boundary patch"
        << abort(FatalError);

    return -1; // Never reached
}

Foam::label Foam::dgGeomMesh::getLocalFaceID
(
    const label faceID,
    const label patchID
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const polyPatch& pp = patches[patchID];

    const label start = pp.start();
    const label size  = pp.size();

    if (faceID < start || faceID >= start + size)
    {
        FatalErrorInFunction
            << "Face ID " << faceID << " is not within patch " << patchID
            << " [start=" << start << ", size=" << size << "]." << nl
            << "→ Valid face IDs for this patch are in range ["
            << start << ", " << start + size - 1 << "]."
            << abort(FatalError);
    }

    return faceID - start;
}

// ************************************************************************* //
