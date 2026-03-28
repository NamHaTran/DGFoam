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

#include "dgGeomMesh.H"
#include "dgBasisField.H"
#include "dgGeomFace.H"
#include "dgRefFace.H"
#include "dgGeomCell.H"

#include "faceList.H"
#include "cellList.H"

#include "polyBoundaryMesh.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "DynamicList.H"

#include "IOobject.H"
#include "IOdictionary.H"
#include "Pstream.H"
#include "error.H"

namespace
{

Foam::label readPOrder(const Foam::fvMesh& mesh)
{
    Foam::IOdictionary dgSchemesDict
    (
        Foam::IOobject
        (
            "dgSchemes",
            mesh.time().system(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        )
    );

    if (!dgSchemesDict.found("dgDiscretization"))
    {
        FatalIOErrorInFunction(dgSchemesDict)
            << "Missing sub-dictionary 'dgDiscretization' in system/dgSchemes."
            << Foam::nl << exit(Foam::FatalIOError);
    }

    const Foam::dictionary& dgDiscretizationDict =
        dgSchemesDict.subDict("dgDiscretization");

    if (!dgDiscretizationDict.found("pOrder"))
    {
        FatalIOErrorInFunction(dgDiscretizationDict)
            << "Missing entry 'pOrder' in system/dgSchemes/dgDiscretization."
            << Foam::nl << exit(Foam::FatalIOError);
    }

    const Foam::label pOrder = dgDiscretizationDict.get<Foam::label>("pOrder");

    if (pOrder < 0)
    {
        FatalIOErrorInFunction(dgDiscretizationDict)
            << "Entry 'pOrder' must be non-negative, but found "
            << pOrder << '.'
            << Foam::nl << exit(Foam::FatalIOError);
    }

    return pOrder;
}

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from fvMesh and read polynomial order from system/dgSchemes
Foam::dgGeomMesh::dgGeomMesh
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    pOrder_(readPOrder(mesh)),
    refFace_(std::make_shared<dgRefFace>(pOrder_)),
    refCellTet_(nullptr),
    refCellHex_(nullptr),
    refCellPrism_(nullptr),
    refCellPyramid_(nullptr)
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
            if (!refCellTet_)
            {
                refCellTet_ = std::make_shared<dgRefCell>(pOrder_, dgCellType::TET);
            }
            refCell = refCellTet_;
        }
        else if (nPoints == 8)
        {
            if (!refCellHex_)
            {
                refCellHex_ = std::make_shared<dgRefCell>(pOrder_, dgCellType::HEX);
            }
            refCell = refCellHex_;
        }
        else if (nPoints == 6)
        {
            if (!refCellPrism_)
            {
                refCellPrism_ = std::make_shared<dgRefCell>(pOrder_, dgCellType::PRISM);
            }
            refCell = refCellPrism_;
        }
        else if (nPoints == 5)
        {
            if (!refCellPyramid_)
            {
                refCellPyramid_ = std::make_shared<dgRefCell>(pOrder_, dgCellType::PYRAMID);
            }
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

        // Calculate Lame parameters at face Gauss points for both owner and neighbour cells
        faces_[faceI]->computeLameParameters();
    }

    // Get boundary faces
    getBoundaryFaces();

    // Assign face connectivity to processor patches if running in parallel
    assignFaceConnectivityToProcPatches();

    // Basis data depend only on geometry/p-order, so cache them once here.
    basisFields_.setSize(mesh_.nCells());
    for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
    {
        basisFields_[cellI] = new dgBasisField(cellI, *this);
    }
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

    forAll(basisFields_, i)
    {
        delete basisFields_[i];
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

        // Skip patches of type 'empty' (used in 2D)
        if (isA<emptyPolyPatch>(pp))
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

        // Skip 'empty' patch type
        if (isA<emptyPolyPatch>(pp))
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

void Foam::dgGeomMesh::assignFaceConnectivityToProcPatches()
{
    // Read dgFacePhysicGauss when running in parallel
    //List<List<vector>> dgFacePhysicGauss;
    //List<List<scalar>> dgFaceJacobian;

    if (Foam::Pstream::parRun())
    {
        IOList<List<vector>> dgFacePhysicGauss
        (
            IOobject
            (
                "dgFacePhysicGauss",
                mesh_.time().constant()/"polyMesh",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        IOList<List<scalar>> dgFaceJacobian
        (
            IOobject
            (
                "dgFaceJacobian",
                mesh_.time().constant()/"polyMesh",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Get polyBoundaryMesh
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Loop over processor patches
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            // Skip patches of type 'empty' (used in 2D)
            if (!isA<processorPolyPatch>(pp))
            {
                continue;
            }

            // Get the starting global face ID of this patch
            const label start = pp.start();

            // Get the number of faces in this patch
            const label size = pp.size();

            // Assign connectivity to each face on this processor patch
            for (label i = 0; i < size; ++i)
            {
                const label faceID = start + i;
                faces_[faceID]->setProcessorPatch(true);
                faces_[faceID]->setPhysicGaussPoints(dgFacePhysicGauss[faceID]);
                faces_[faceID]->setLameParams(dgFaceJacobian[faceID]);
                faces_[faceID]->calcGaussPointsOnOwnerSide();

                // Processor faces inherit a canonical physical Gauss-point
                // ordering from prepareDGDecomposedMesh. Rebuild owner-side
                // basis/derivative data on that updated ordering so flux
                // reconstruction and communicated traces use the same face
                // quadrature layout on every rank.
                faces_[faceID]->computeBasisAndDerivatives();
            }
        }
    }
}


const Foam::dgBasisField& Foam::dgGeomMesh::basisField(const label cellID) const
{
    if (cellID < 0 || cellID >= basisFields_.size() || !basisFields_[cellID])
    {
        FatalErrorInFunction
            << "Invalid cached dgBasisField access for cell " << cellID
            << abort(FatalError);
    }

    return *basisFields_[cellID];
}
// ************************************************************************* //
