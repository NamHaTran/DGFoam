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

#include "dgGeomFace.H"
#include "pointField.H"
#include "fvMesh.H"
#include "vector2D.H"
#include "SortList.H"
#include "face.H"
#include "pointField.H"
#include "mathematicalConstants.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Default constructor
Foam::dgGeomFace::dgGeomFace()
:
    faceID_(-1),
    mesh_(*static_cast<fvMesh*>(nullptr)), // placeholder invalid ref
    refFace_(nullptr),
    patchID_(-1),
    isBoundary_(false),
    isProcessorPatch_(false)
{}

// Construct from face and mesh
Foam::dgGeomFace::dgGeomFace
(
    label faceID,
    const fvMesh& mesh,
    const std::shared_ptr<dgRefFace>& refFace
)
:
    faceID_(faceID),
    mesh_(mesh),
    refFace_(refFace),
    patchID_(-1),
    isBoundary_(false),
    isProcessorPatch_(false)
{
    // Get face from faceID_
    const face& f = mesh_.faces()[faceID_];

    // Get number of face node
    const label nPoints = f.size();

    // Get face points in original order
    globalPoints_.setSize(nPoints);

    // Get 3D coordinates of face points
    for (label i = 0; i < nPoints; ++i)
    {
        globalPoints_[i] = mesh_.points()[f[i]];
    }

    // Get boundary information (if any)
    patchID_ = Foam::findOwnerPatch(mesh_, faceID_);
    isBoundary_ = (patchID_ != -1);

    processFlatAndSortedPoints();

    // Detect face type
    if (nPoints == 3)
    {
        type_ = dgFaceType::TRI;
    }
    else if (nPoints == 4)
    {
        type_ = dgFaceType::QUAD;
    }
    else
    {
        FatalErrorInFunction
            << "Face must have at least 3 points." << abort(FatalError);
    }
}

// Copy constructor
Foam::dgGeomFace::dgGeomFace(const dgGeomFace& other)
:
    faceID_(other.faceID_),
    type_(other.type_),
    mesh_(other.mesh_),
    refFace_(other.refFace_),
    globalPoints_(other.globalPoints_),
    ownerPos_(other.ownerPos_),
    neighborPos_(other.neighborPos_),
    ownerCellType_(other.ownerCellType_),
    neighborCellType_(other.neighborCellType_),
    flattenedPoints_(other.flattenedPoints_),
    connectivity_(other.connectivity_),
    ownerBasisData_(other.ownerBasisData_),
    neighborBasisData_(other.neighborBasisData_),
    ownerJ2D_(other.ownerJ2D_),
    neighborJ2D_(other.neighborJ2D_),
    patchID_(other.patchID_),
    isBoundary_(other.isBoundary_),
    isProcessorPatch_(other.isProcessorPatch_)
{

}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return number of vertices in the face
Foam::label Foam::dgGeomFace::size() const
{
    return mesh_.faces()[faceID_].size();
}


// Return base face (index list into mesh.points())
const Foam::face& Foam::dgGeomFace::baseFace() const
{
    return mesh_.faces()[faceID_];
}


// Return coordinate of the i-th point on the face
const Foam::point& Foam::dgGeomFace::getPoint(label i) const
{
    const face& f = mesh_.faces()[faceID_];
    return mesh_.points()[f[i]];
}


// Return the face center (using foam::face::centre)
Foam::point Foam::dgGeomFace::centre() const
{
    //return mesh_.Cf()[faceID_];
    return mesh_.faceCentres()[faceID_];
}


Foam::vector Foam::dgGeomFace::areaNormal() const
{
    //return mesh_.Sf()[faceID_];
    return mesh_.faceAreas()[faceID_];
}


// Return face normal (assuming planar face)
Foam::vector Foam::dgGeomFace::normal() const
{
    return mesh_.faceAreas()[faceID_] / mag(mesh_.faceAreas()[faceID_]);
}


// Return face area
Foam::scalar Foam::dgGeomFace::area() const
{
    return mag(mesh_.faceAreas()[faceID_]);
}


// Copy assignment operator
Foam::dgGeomFace& Foam::dgGeomFace::operator=(const dgGeomFace& other)
{
    if (this != &other)
    {
        faceID_           = other.faceID_;
        refFace_          = other.refFace_;
        globalPoints_     = other.globalPoints_;
        ownerPos_         = other.ownerPos_;
        neighborPos_      = other.neighborPos_;
        ownerCellType_    = other.ownerCellType_;
        neighborCellType_ = other.neighborCellType_;
        // sortedPointLabels_ = other.sortedPointLabels_;
        flattenedPoints_  = other.flattenedPoints_;
        type_             = other.type_;
        connectivity_     = other.connectivity_;
        ownerBasisData_   = other.ownerBasisData_;
        neighborBasisData_ = other.neighborBasisData_;
        ownerJ2D_         = other.ownerJ2D_;
        neighborJ2D_      = other.neighborJ2D_;
        patchID_          = other.patchID_;
        isBoundary_       = other.isBoundary_;
        isProcessorPatch_ = other.isProcessorPatch_;
    }
    return *this;
}

void Foam::dgGeomFace::printDebugInfo() const
{
    const face& f = mesh_.faces()[faceID_];
    const pointField& points = mesh_.points();

    Info << "[dgGeomFace::printDebugInfo()] Face ID: " << faceID_ << nl;
    Info << "  Number of points: " << f.size() << " (original order):" << nl;

    forAll(f, i)
    {
        label ptID = f[i];
        const point& pt = points[ptID];
        Info << "    Point " << i << " (ID " << ptID << "): " << pt << nl;
    }

    Info << "  Type                 : " << type_ << nl;
    Info << "  Owner position       : " << ownerPos_ << nl;
    Info << "  Neighbor position    : " << neighborPos_ << nl;
    Info << "  Owner cell type      : " << ownerCellType_ << nl;
    Info << "  Neighbor cell type   : " << neighborCellType_ << nl;

    Info << "  Global points (CCW from owner):" << nl;
    forAll(globalPoints_, i)
    {
        Info << "    " << globalPoints_[i] << nl;
    }

    Info << "  Flattened points:" << nl;
    forAll(flattenedPoints_, i)
    {
        Info << "    " << flattenedPoints_[i] << nl;
    }

    Info << "  Gauss connectivity (owner → neighbor):" << nl;
    forAll(connectivity_, i)
    {
        Info << "    owner[" << i << "] → neighbor[" << connectivity_[i] << "]" << nl;
    }

    Info << "  Basis functions at Gauss points (owner side):" << nl;
    {
        for (label gp = 0; gp < ownerBasisData_.basis.size(); ++gp)
        {
            Info << "    Gauss point " << gp << ":" << nl;
            Info << "      basis        : " << ownerBasisData_.basis[gp] << nl;
            Info << "      dBasis/dEta1 : " << ownerBasisData_.dBasis_dEta1[gp] << nl;
            Info << "      dBasis/dEta2 : " << ownerBasisData_.dBasis_dEta2[gp] << nl;
            Info << "      dBasis/dEta3 : " << ownerBasisData_.dBasis_dEta3[gp] << nl;
        }
    }

    if (neighborCellType_ != dgCellType::NONE)
    {
        Info << "  Basis functions at Gauss points (neighbor side):" << nl;
        for (label gp = 0; gp < neighborBasisData_.basis.size(); ++gp)
        {
            Info << "    Gauss point " << gp << ":" << nl;
            Info << "      basis        : " << neighborBasisData_.basis[gp] << nl;
            Info << "      dBasis/dEta1 : " << neighborBasisData_.dBasis_dEta1[gp] << nl;
            Info << "      dBasis/dEta2 : " << neighborBasisData_.dBasis_dEta2[gp] << nl;
            Info << "      dBasis/dEta3 : " << neighborBasisData_.dBasis_dEta3[gp] << nl;
        }
    }
}

// Flatten face points to 2D using face normal
void Foam::dgGeomFace::flattenFace()
{
    const label nPoints = globalPoints_.size();
    flattenedPoints_.setSize(nPoints);

    // Compute centroid of face
    vector centroid = vector::zero;
    forAll(globalPoints_, i)
    {
        centroid += globalPoints_[i];
    }
    centroid /= nPoints;

    // Compute face normal
    vector n = normal();
    n /= mag(n) + VSMALL;

    // Create orthonormal basis (u, v) on the face
    vector u;
    if (mag(n ^ vector(1, 0, 0)) > VSMALL)
    {
        u = n ^ vector(1, 0, 0);
    }
    else
    {
        u = n ^ vector(0, 1, 0);
    }
    u /= mag(u) + VSMALL;

    vector v = n ^ u;
    v /= mag(v) + VSMALL;

    // Project each point onto (u,v) plane centered at centroid
    for (label i = 0; i < nPoints; ++i)
    {
        vector r = globalPoints_[i] - centroid;
        scalar x = (r & u);
        scalar y = (r & v);
        flattenedPoints_[i] = vector2D(x, y);
    }
}

// Sort flattened face points in counter-clockwise order
/*
void Foam::dgGeomFace::sortPointsCCW()
{
    const label nPoints = flattenedPoints_.size();
    if (nPoints < 3)
    {
        FatalErrorInFunction
            << "Face must have at least 3 points for CCW sorting." << abort(FatalError);
    }

    // Compute centroid
    vector2D centroid(0, 0);
    for (const auto& pt : flattenedPoints_)
    {
        centroid += pt;
    }
    centroid /= scalar(nPoints);

    // Compute angle between point and centroid
    List<scalar> angles(nPoints);
    for (label i = 0; i < nPoints; ++i)
    {
        const scalar dx = flattenedPoints_[i].x() - centroid.x();
        const scalar dy = flattenedPoints_[i].y() - centroid.y();
        angles[i] = atan2(dy, dx);
    }

    // Sort by angle using SortList
    SortList<scalar> sortedAngles(angles);

    // Reorder point labels and flattened coords
    const face& f = mesh_.faces()[faceID_];
    sortedPointLabels_.setSize(nPoints);
    List<vector2D> sortedFlat(nPoints);

    for (label i = 0; i < nPoints; ++i)
    {
        label oldI = sortedAngles.indices()[i];
        sortedPointLabels_[i] = f[oldI];         // global point ID
        sortedFlat[i]         = flattenedPoints_[oldI];
    }

    flattenedPoints_ = sortedFlat;
}
*/

// Perform flattening, sorting and store data
void Foam::dgGeomFace::processFlatAndSortedPoints()
{
    flattenFace();
    // sortPointsCCW();
}

// Mapping from reference to physical space for quad or tri face
void Foam::dgGeomFace::mappingFromRefToReal
(
    const dgFaceType type,
    const List<vector>& gaussPoints,
    const List<vector2D>& faceVertices,
    List<vector2D>& physicGaussP
)
{
    physicGaussP.setSize(gaussPoints.size());

    forAll(gaussPoints, i)
    {
        scalar eta1 = gaussPoints[i].x();
        scalar eta2 = gaussPoints[i].y();
        scalar eta3 = gaussPoints[i].z();

        scalar eta1_flat = 0;
        scalar eta2_flat = 0;

        // Remove coordinate equal to 1 or -1 to get 2D projection
        if (mag(eta1) == 1.0)
        {
            eta1_flat = eta2;
            eta2_flat = eta3;
        }
        else if (mag(eta2) == 1.0)
        {
            eta1_flat = eta1;
            eta2_flat = eta3;
        }
        else if (mag(eta3) == 1.0)
        {
            eta1_flat = eta1;
            eta2_flat = eta2;
        }

        if (type == dgFaceType::QUAD)
        {
            scalar C1 = 0.25 * (1.0 - eta1_flat) * (1.0 - eta2_flat);
            scalar C2 = 0.25 * (1.0 + eta1_flat) * (1.0 - eta2_flat);
            scalar C3 = 0.25 * (1.0 + eta1_flat) * (1.0 + eta2_flat);
            scalar C4 = 0.25 * (1.0 - eta1_flat) * (1.0 + eta2_flat);

            const vector2D& v1 = faceVertices[0];
            const vector2D& v2 = faceVertices[1];
            const vector2D& v3 = faceVertices[2];
            const vector2D& v4 = faceVertices[3];

            scalar x = C1 * v1.x() + C2 * v2.x() + C3 * v3.x() + C4 * v4.x();
            scalar y = C1 * v1.y() + C2 * v2.y() + C3 * v3.y() + C4 * v4.y();

            physicGaussP[i] = vector2D(x, y);
        }
        else if (type == dgFaceType::TRI)
        {
            scalar C1 = 0.25 * (1.0 - eta1_flat) * (1.0 - eta2_flat);
            scalar C2 = 0.25 * (1.0 + eta1_flat) * (1.0 - eta2_flat);
            scalar C3 = 0.5  * (1.0 + eta2_flat);

            const vector2D& v1 = faceVertices[0];
            const vector2D& v2 = faceVertices[1];
            const vector2D& v3 = faceVertices[2];

            scalar x = C1 * v1.x() + C2 * v2.x() + C3 * v3.x();
            scalar y = C1 * v1.y() + C2 * v2.y() + C3 * v3.y();

            physicGaussP[i] = vector2D(x, y);
        }
        else
        {
            FatalErrorInFunction << "Unsupported face type." << abort(FatalError);
        }
    }
}

// Compute connectivity between Gauss points from owner and neighbor
void Foam::dgGeomFace::findGaussConnectivity()
{
    const List<vector>& gaussOwner = gaussPointsOwner();
    const label nGauss = gaussOwner.size();
    connectivity_.setSize(nGauss);

    if (neighborPos_ == dgFacePosition::NONE)
    {
        for (label i = 0; i < nGauss; ++i)
        {
            connectivity_[i] = i;
        }
        return;
    }

    // Step 1: generate flippedPoints from flattenedPoints_
    const label nPoints = flattenedPoints_.size();
    List<vector2D> flippedPoints(nPoints);
    for (label i = 0; i < nPoints; ++i)
    {
        flippedPoints[i] = flattenedPoints_[nPoints - 1 - i];
    }

    // Step 2: get Gauss points on neighbor side
    const List<vector>& gaussNeighbor = gaussPointsNeighbor();

    // Step 3: map to physical space
    List<vector2D> physicGaussOwner;
    List<vector2D> physicGaussNeighbor;

    mappingFromRefToReal(type_, gaussOwner, flattenedPoints_, physicGaussOwner);
    mappingFromRefToReal(type_, gaussNeighbor, flippedPoints, physicGaussNeighbor);

    // Step 4: compare and find connectivity
    const scalar tol = 1e-10;

    for (label i = 0; i < nGauss; ++i)
    {
        const vector2D& pO = physicGaussOwner[i];

        for (label j = 0; j < nGauss; ++j)
        {
            const vector2D& pN = physicGaussNeighbor[j];

            if
            (
                mag(pO.x() - pN.x()) < tol &&
                mag(pO.y() - pN.y()) < tol
            )
            {
                connectivity_[i] = j;
                break;
            }
        }
    }
}

void Foam::dgGeomFace::computeBasisAndDerivatives()
{
    // Calculate basis functions and derivatives for both owner and neighbor
    ownerBasisData_ = refFace_->computeBasisAndDerivatives(ownerCellType_, ownerPos_);

    // Only compute neighbor if it exists
    if (neighborPos_ != dgFacePosition::NONE)
    {
        neighborBasisData_ = refFace_->computeBasisAndDerivatives(neighborCellType_, neighborPos_);
    }
}

void Foam::dgGeomFace::computeOwnerLameParameters
(
    const List<vector>& cellVertices
)
{
    // Calculate Jacobian 2D at Gauss points
    const List<vector>& ownerGaussPts = refFace_->gaussPoints(ownerPos_);
    const label nGauss = ownerGaussPts.size();
    ownerJ2D_.setSize(nGauss);

    for (label gp = 0; gp < nGauss; ++gp)
    {
        ownerJ2D_[gp] = calcLameParam(ownerCellType_, ownerPos_, ownerGaussPts[gp], cellVertices);
    }
}

void Foam::dgGeomFace::computeNeighborLameParameters
(
    const List<vector>& cellVertices
)
{
    // Calculate Jacobian 2D at Gauss points
    const List<vector>& neighborGaussPts = refFace_->gaussPoints(neighborPos_);
    const label nGauss = neighborGaussPts.size();
    neighborJ2D_.setSize(nGauss);

    for (label gp = 0; gp < nGauss; ++gp)
    {
        neighborJ2D_[gp] = calcLameParam(neighborCellType_, neighborPos_, neighborGaussPts[gp], cellVertices);
    }
}

// ************************************************************************* //
