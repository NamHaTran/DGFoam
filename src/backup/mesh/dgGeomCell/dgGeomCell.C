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

#include "dgGeomCell.H"
#include "faceList.H"
#include "cellList.H"
#include "cellShape.H"
#include "HashSet.H"
#include "Jacobian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dgGeomCell::dgGeomCell
(
    const label cellID,
    const fvMesh& mesh,
    const std::shared_ptr<dgRefCell>& refCell
)
:
    cellID_(cellID),
    mesh_(mesh),
    refCell_(refCell)
{
    // Access cell shape
    const cellShape& shape = mesh_.cellShapes()[cellID_];

    // Get cell vertices in standard order
    cellPoints_ = shape.points(mesh_.points());

    // Set cell type
    const label nPoints = cellPoints_.size();

    if (nPoints == 4)
    {
        type_ = dgCellType::TET;
    }
    else if (nPoints == 8)
    {
        type_ = dgCellType::HEX;
    }
    else if (nPoints == 6)
    {
        type_ = dgCellType::PRISM;
    }
    else if (nPoints == 5)
    {
        type_ = dgCellType::PYRAMID;
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported cell type with " << nPoints << " points." << abort(FatalError);
    }

    // Get face labels in model order (global face IDs)
    faceLabels_ = shape.meshFaces(mesh_.faces(), mesh_.cells()[cellID_]);

    // Calculate internal Jacobians at Gauss points
    const List<vector>& gaussPts = refCell_->gaussPoints();
    const label nGauss = gaussPts.size();
    J3D_.setSize(nGauss);

    for (label gp = 0; gp < nGauss; ++gp)
    {
        J3D_[gp] = calcJacobianDetAtInteriorGaussPt(type_, gaussPts[gp], cellPoints_);
    }

    // Precompute mass matrix
    const List<List<scalar>>& basis = refCell_->basis();

    // Gauss weights
    const List<scalar>& w = refCell_->weights();

    const label nDof = basis[0].size();

    // Resize and initialize mass matrix
    massMatrix_.resize(nDof);

    // Assemble mass matrix directly at Gauss points
    for (label i = 0; i < nDof; ++i)
    {
        for (label j = 0; j < nDof; ++j)
        {
            scalar mij = 0.0;

            for (label gp = 0; gp < nGauss; ++gp)
            {
                mij +=
                    basis[gp][i]
                  * basis[gp][j]
                  * J3D_[gp]
                  * w[gp];
            }

            massMatrix_[i][j] = mij;
        }
    }
}

Foam::dgGeomCell::~dgGeomCell()
{
    //forAll(faces_, i)
    //{
    //    delete faces_[i];
    //    faces_[i] = nullptr;
    //}
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return the centroid of the cell.
// Note:
// - mesh_.C() is a volVectorField (GeometricField)
// - We must use .internalField()[cellID_] to access internal value
Foam::point Foam::dgGeomCell::centre() const
{
    return mesh_.C().internalField()[cellID_];
}

// Return the volume of the cell.
// Note:
// - mesh_.V() is a DimensionedField<scalar, volMesh>
// - It supports direct access via [cellID_] (no .internalField() needed)
Foam::scalar Foam::dgGeomCell::volume() const
{
    return mesh_.V()[cellID_];
}


void Foam::dgGeomCell::printDebugInfo() const
{
    Info << "\n========== Cell " << cellID_ << " ==========" << nl;
    Info << " - Centre: " << centre() << nl;
    Info << " - Volume: " << volume() << nl;

    Info << " - Number of faces: " << nFaces() << nl;
    Info << " - Number of points: " << nPoints() << nl;
    forAll(cellPoints_, i)
    {
        Info << "     Point local ID " << i << ": " << cellPoints_[i] << nl;
    }

    const List<vector>& gaussPts = gaussPoints();
    const List<scalar>& weights  = this->weights();

    const List<List<scalar>>& b   = basis();
    const List<List<scalar>>& db1 = dBasis_dEta1();
    const List<List<scalar>>& db2 = dBasis_dEta2();
    const List<List<scalar>>& db3 = dBasis_dEta3();

    const label nGauss = gaussPts.size();
    const label nBasis = b[0].size();

    Info << " - Number of Gauss points: " << nGauss << nl;

    for (label gp = 0; gp < nGauss; ++gp)
    {
        Info << "   Gauss Point [" << gp << "]" << nl;
        Info << "     eta       : " << gaussPts[gp] << nl;
        Info << "     weight    : " << weights[gp] << nl;

        Info << "     basis     :";
        for (label k = 0; k < nBasis; ++k)
        {
            Info << " " << b[gp][k];
        }
        Info << nl;

        Info << "     dBasis/dη1:";
        for (label k = 0; k < nBasis; ++k)
        {
            Info << " " << db1[gp][k];
        }
        Info << nl;

        Info << "     dBasis/dη2:";
        for (label k = 0; k < nBasis; ++k)
        {
            Info << " " << db2[gp][k];
        }
        Info << nl;

        Info << "     dBasis/dη3:";
        for (label k = 0; k < nBasis; ++k)
        {
            Info << " " << db3[gp][k];
        }
        Info << nl;

        Info << "     Jacobian det:";
        Info << " " << J3D_[gp];
        Info << nl;
    }

    Info << " - Mass matrix:" << massMatrix_ << nl;

    Info << "==========================================" << endl;
}

//- Update faces information based on the cell's geometry and mesh
void Foam::dgGeomCell::updateFaceInfo
(
    List<dgGeomFace*>& faces
)
{
    neighborCellLabels_.resize(faceLabels_.size());

    forAll(faceLabels_, localID)
    {
        const label faceI = faceLabels_[localID];

        // Sanity check
        if (faceI < 0 || faceI >= faces.size())
        {
            FatalErrorInFunction
                << "Invalid face index " << faceI << " for cell " << cellID_ << nl
                << abort(FatalError);
        }

        dgGeomFace* facePtr = faces[faceI];

        // Convert local face ID to reference position
        dgFacePosition pos;

        switch (type_)
        {
            case dgCellType::HEX:
                pos = convertIdToPositionOnHex(localID);
                break;

            case dgCellType::PRISM:
                pos = mapFacePositionFromPrism(convertIdToPositionOnPrism(localID));
                break;

            case dgCellType::TET:
                pos = mapFacePositionFromTet(convertIdToPositionOnTet(localID));
                break;

            case dgCellType::PYRAMID:
                pos = mapFacePositionFromPyramid(convertIdToPositionOnPyramid(localID));
                break;

            default:
                FatalErrorInFunction
                    << "Unsupported cell type: " << type_ << nl
                    << abort(FatalError);
        }

        // Determine ownership and boundary condition
        const label owner = mesh_.faceOwner()[faceI];
        const bool isBoundaryFace = (faceI >= mesh_.nInternalFaces());

        if (owner == cellID_)
        {
            facePtr->setOwnerPos(pos);
            facePtr->setOwnerCellType(type_);
            
            // Calculate Lame parameters at face Gauss points
            facePtr->computeOwnerLameParameters(cellPoints_);

            if (isBoundaryFace)
            {
                facePtr->setNeighborPos(dgFacePosition::NONE);
                facePtr->setNeighborCellType(dgCellType::NONE);
                neighborCellLabels_[localID] = -1;
            }
            else
            {
                const label neighbor = mesh_.faceNeighbour()[faceI];
                neighborCellLabels_[localID] = neighbor;
            }
        }
        else
        {
            facePtr->setNeighborPos(pos);
            facePtr->setNeighborCellType(type_);
            neighborCellLabels_[localID] = owner;

            // Calculate Lame parameters at face Gauss points
            facePtr->computeNeighborLameParameters(cellPoints_);
        }
    }
}

// ************************************************************************* //

