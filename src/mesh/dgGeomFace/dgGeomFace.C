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

#include "dgGeomFace.H"
#include "pointField.H"
#include "fvMesh.H"
#include "vector2D.H"
#include "SortList.H"
#include "face.H"
#include "pointField.H"
#include "mathematicalConstants.H"
#include "dgAffineMapping.H"
#include "refCoordTransforms.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Default constructor
Foam::dgGeomFace::dgGeomFace()
:
    faceID_(-1),
    type_(dgFaceType::INVALID),
    mesh_(*static_cast<fvMesh*>(nullptr)), // placeholder invalid ref
    refFace_(nullptr),
    nGauss_(0),
    ownerId_(-1),
    neighborId_(-1),
    ownerPos_(dgFacePosition::NONE),
    neighborPos_(dgFacePosition::NONE),
    ownerCellType_(dgCellType::NONE),
    neighborCellType_(dgCellType::NONE),
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
    type_(dgFaceType::INVALID),
    mesh_(mesh),
    refFace_(refFace),
    nGauss_(0),
    ownerId_(-1),
    neighborId_(-1),
    ownerPos_(dgFacePosition::NONE),
    neighborPos_(dgFacePosition::NONE),
    ownerCellType_(dgCellType::NONE),
    neighborCellType_(dgCellType::NONE),
    patchID_(-1),
    isBoundary_(false),
    isProcessorPatch_(false)
{
    // Get Number of Gauss point
    nGauss_ = refFace_->nGauss();

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
    physGaussPoints_(other.physGaussPoints_),
    gaussPointsOwner_(other.gaussPointsOwner_),
    gaussPointsNeighbor_(other.gaussPointsNeighbor_),
    ownerBasisData_(other.ownerBasisData_),
    neighborBasisData_(other.neighborBasisData_),
    J2D_(other.J2D_),
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
        physGaussPoints_  = other.physGaussPoints_;
        gaussPointsOwner_  = other.gaussPointsOwner_;
        gaussPointsNeighbor_ = other.gaussPointsNeighbor_;
        type_             = other.type_;
        ownerBasisData_   = other.ownerBasisData_;
        neighborBasisData_ = other.neighborBasisData_;
        J2D_              = other.J2D_;
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

    Pout << "[dgGeomFace::printDebugInfo()] Face ID: " << faceID_ << nl;
    Pout << "  Number of points: " << f.size() << " (original order):" << nl;

    if (this->isBoundary())
    {
        Pout << "  This face is a boundary face on patch " << patchID_ << nl;
    }
    else
    {
        Pout << "  This face is an internal face." << nl;
    }

    forAll(f, i)
    {
        label ptID = f[i];
        const point& pt = points[ptID];
        Pout << "    Point " << i << " (ID " << ptID << "): " << pt << nl;
    }

    Pout << "  Type                 : " << type_ << nl;
    Pout << "  Owner position       : " << ownerPos_ << nl;
    Pout << "  Neighbor position    : " << neighborPos_ << nl;
    Pout << "  Owner cell type      : " << ownerCellType_ << nl;
    Pout << "  Neighbor cell type   : " << neighborCellType_ << nl;

    Pout << "  Global points (CCW from owner):" << nl;
    forAll(globalPoints_, i)
    {
        Pout << "    " << globalPoints_[i] << nl;
    }

    Pout << "  Flattened points:" << nl;
    forAll(flattenedPoints_, i)
    {
        Pout << "    " << flattenedPoints_[i] << nl;
    }

    Pout << "  Physical Gauss points:" << nl;
    forAll(physGaussPoints_, i)
    {        
        Pout << "    " << physGaussPoints_[i] << nl;
    }

    Pout << "  2D Gauss points on reference face:" << nl;
    forAll(refFace_->gaussPoints(), i)
    {        
        Pout << "    " << refFace_->gaussPoints()[i] << nl;
    }

    Pout << "  Gauss points on owner reference face:" << nl;
    forAll(gaussPointsOwner_, i)
    {        
        Pout << "    " << gaussPointsOwner_[i] << nl;
    }

    Pout << "  Gauss points on neighbor reference face:" << nl;
    forAll(gaussPointsNeighbor_, i)
    {        
        Pout << "    " << gaussPointsNeighbor_[i] << nl;
    }

    Pout << "  Basis functions at Gauss points (owner side):" << nl;
    {
        for (label gp = 0; gp < ownerBasisData_.basis.size(); ++gp)
        {
            Pout << "    Gauss point " << gp << ":" << nl;
            Pout << "      basis        : " << ownerBasisData_.basis[gp] << nl;
            Pout << "      dBasis/dEta1 : " << ownerBasisData_.dBasis_dEta1[gp] << nl;
            Pout << "      dBasis/dEta2 : " << ownerBasisData_.dBasis_dEta2[gp] << nl;
            Pout << "      dBasis/dEta3 : " << ownerBasisData_.dBasis_dEta3[gp] << nl;
        }
    }

    if (neighborCellType_ != dgCellType::NONE)
    {
        Pout << "  Basis functions at Gauss points (neighbor side):" << nl;
        for (label gp = 0; gp < neighborBasisData_.basis.size(); ++gp)
        {
            Pout << "    Gauss point " << gp << ":" << nl;
            Pout << "      basis        : " << neighborBasisData_.basis[gp] << nl;
            Pout << "      dBasis/dEta1 : " << neighborBasisData_.dBasis_dEta1[gp] << nl;
            Pout << "      dBasis/dEta2 : " << neighborBasisData_.dBasis_dEta2[gp] << nl;
            Pout << "      dBasis/dEta3 : " << neighborBasisData_.dBasis_dEta3[gp] << nl;
        }
    }

    Pout << "  Lame parameters at Gauss points:" << nl;
    forAll(J2D_, i)
    {        
        Pout << "    Gauss point " << i << ": J2D = " << J2D_[i] << nl;
    }
}

void Foam::dgGeomFace::printCellsInfo() const
{
    if (neighborId_ != -1)
    {
        // Cells information
        // Access cell shape
        const cellShape& ownerShape = mesh_.cellShapes()[ownerId_];

        // Get cell vertices in standard order
        pointField ownerCellPoints = ownerShape.points(mesh_.points());

        const cellShape& neighborShape = mesh_.cellShapes()[neighborId_];
        pointField neighborCellPoints = neighborShape.points(mesh_.points());

        labelList ownerFaceLabels_ = ownerShape.meshFaces(mesh_.faces(), mesh_.cells()[ownerId_]);
        labelList neighborFaceLabels_ = neighborShape.meshFaces(mesh_.faces(), mesh_.cells()[neighborId_]);

        Info << "Owner cell ID: " << ownerId_ << ", type: " << ownerCellType_ << nl;
        Info << "Neighbor cell ID: " << neighborId_ << ", type: " << neighborCellType_ << nl;
        Info << "Owner cell points:" << ownerCellPoints << nl;
        Info << "Neighbor cell points:" << neighborCellPoints << nl;
        Info << "Owner cell face labels:" << ownerFaceLabels_ << nl;
        Info << "Neighbor cell face labels:" << neighborFaceLabels_ << nl;

        Info << "  Check mapping consistency:" << nl;
        forAll(gaussPointsOwner_, i)
        {
            vector physicPt = Foam::mapEtaToX(gaussPointsOwner_[i], ownerCellPoints, ownerCellType_);
            Info << "    Owner Gauss point " << i << ": physPt = " << physicPt << nl;
        }

        if (neighborId_ != -1)
        {
            forAll(gaussPointsNeighbor_, i)
            {
                vector physicPt = Foam::mapEtaToX(gaussPointsNeighbor_[i], neighborCellPoints, neighborCellType_);
                Info << "    Neighbor Gauss point " << i << ": physPt = " << physicPt << nl;
            }
        }
    }
}

// Flatten the face points to 2D coordinates in the face plane
void Foam::dgGeomFace::flattenFace
(
    const List<vector>& globalPoints,   // Input 3D coordinates
    List<vector2D>& flattenedPoints     // Output 2D coordinates
) const
{
    // ---------------------------------------------------------------------
    // Number of points
    // ---------------------------------------------------------------------

    const label nPoints = globalPoints.size();

    flattenedPoints.setSize(nPoints);


    // ---------------------------------------------------------------------
    // Compute centroid of the face
    // ---------------------------------------------------------------------

    vector centroid = vector::zero;

    forAll(globalPoints, i)
    {
        centroid += globalPoints[i];
    }

    centroid /= nPoints;


    // ---------------------------------------------------------------------
    // Compute normalized face normal
    // ---------------------------------------------------------------------

    vector n = normal();

    // Normalize and avoid division by zero
    n /= mag(n) + VSMALL;


    // ---------------------------------------------------------------------
    // Build orthonormal basis (u, v) on the face
    // ---------------------------------------------------------------------

    vector u(Zero);

    // Choose a reference vector not parallel to n
    if (mag(n ^ vector(1, 0, 0)) > VSMALL)
    {
        u = n ^ vector(1, 0, 0);
    }
    else
    {
        u = n ^ vector(0, 1, 0);
    }

    u /= mag(u) + VSMALL;

    // v completes right-handed orthonormal system
    vector v = n ^ u;
    v /= mag(v) + VSMALL;


    // ---------------------------------------------------------------------
    // Project each point onto (u,v) plane
    // ---------------------------------------------------------------------

    for (label i = 0; i < nPoints; ++i)
    {
        const vector r = globalPoints[i] - centroid;

        const scalar x = (r & u);   // Projection onto u
        const scalar y = (r & v);   // Projection onto v

        flattenedPoints[i] = vector2D(x, y);
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
    flattenedPoints_.setSize(globalPoints_.size());
    flattenFace(globalPoints_, flattenedPoints_);
}

// Compute connectivity between Gauss points from owner and neighbor cells
void Foam::dgGeomFace::findGaussConnectivity()
{
    // Get Gauss points 2D on the reference quad face
    const List<vector2D>& gauss2D = refFace_->gaussPoints();

    // Map to physical space
    physGaussPoints_.setSize(gauss2D.size());

    for (label i = 0; i < nGauss_; ++i)
    {
        if (type_ == dgFaceType::TRI)
        {
            // Calculate barycentric coordinates of the Gauss point in the reference face
            vector2D rs = Foam::squareToTriangle(gauss2D[i]);

            // Map to physical space using the actual face vertices
            physGaussPoints_[i] = Foam::mapToPhysicalTriangle
            (
                rs,
                globalPoints_
            );
        }
        else if (type_ == dgFaceType::QUAD)
        {
            // Use the original 2D Gauss points for quadrilateral faces
            physGaussPoints_[i] = Foam::mapToPhysicalQuad
            (
                gauss2D[i],
                globalPoints_
            );
        }
    }

    // Map Gauss points to owner and neighbor reference space
    gaussPointsOwner_.setSize(nGauss_);
    gaussPointsNeighbor_.setSize(nGauss_);

    // Access cell shape
    const cellShape& ownerShape = mesh_.cellShapes()[ownerId_];
    pointField ownerCellPoints = ownerShape.points(mesh_.points());

    for (label i = 0; i < nGauss_; ++i)
    {
        gaussPointsOwner_[i] = Foam::mapXToEta(physGaussPoints_[i], ownerCellPoints, ownerCellType_);
        if (neighborPos_ != dgFacePosition::NONE)
        {
            const cellShape& neighborShape = mesh_.cellShapes()[neighborId_];
            pointField neighborCellPoints = neighborShape.points(mesh_.points());
            gaussPointsNeighbor_[i] = Foam::mapXToEta(physGaussPoints_[i], neighborCellPoints, neighborCellType_);
        }
        else // Boundary face, no neighbor
        {
            gaussPointsNeighbor_[i] = vector::zero; // or some invalid value
        }
    }
}

void Foam::dgGeomFace::calcGaussPointsOnOwnerSide()
{
    // Map Gauss points to owner and neighbor reference space
    gaussPointsOwner_.setSize(nGauss_);

    // Access cell shape
    const cellShape& ownerShape = mesh_.cellShapes()[ownerId_];
    pointField ownerCellPoints = ownerShape.points(mesh_.points());

    for (label i = 0; i < nGauss_; ++i)
    {
        gaussPointsOwner_[i] = Foam::mapXToEta(physGaussPoints_[i], ownerCellPoints, ownerCellType_);
    }
}

void Foam::dgGeomFace::computeBasisAndDerivatives()
{
    // Calculate basis functions and derivatives for both owner and neighbor
    ownerBasisData_ = refFace_->computeBasisAndDerivatives(gaussPointsOwner_, ownerCellType_);

    // Only compute neighbor if it exists
    if (neighborPos_ != dgFacePosition::NONE)
    {
        neighborBasisData_ = refFace_->computeBasisAndDerivatives(gaussPointsNeighbor_, neighborCellType_);
    }
}

void Foam::dgGeomFace::computeLameParameters()
{
    const List<vector2D>& gauss2D = refFace_->gaussPoints();

    J2D_.setSize(nGauss_);

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        const scalar eta1 = gauss2D[gp].x();
        const scalar eta2 = gauss2D[gp].y();

        if (type_ == dgFaceType::TRI)
        {
            const scalar a = 0.5*(1.0 + eta1);
            const scalar b = 0.5*(1.0 + eta2);

            const vector e01 = globalPoints_[1] - globalPoints_[0];
            const vector e02 = globalPoints_[2] - globalPoints_[0];

            const vector dX_dEta1 = 0.5*(e01 - b*e02);
            const vector dX_dEta2 = 0.5*(1.0 - a)*e02;

            J2D_[gp] = mag(dX_dEta1 ^ dX_dEta2);
        }
        else if (type_ == dgFaceType::QUAD)
        {
            const vector dX_dEta1 =
                -0.25*(1.0 - eta2)*globalPoints_[0]
              +  0.25*(1.0 - eta2)*globalPoints_[1]
              +  0.25*(1.0 + eta2)*globalPoints_[2]
              -  0.25*(1.0 + eta2)*globalPoints_[3];

            const vector dX_dEta2 =
                -0.25*(1.0 - eta1)*globalPoints_[0]
              -  0.25*(1.0 + eta1)*globalPoints_[1]
              +  0.25*(1.0 + eta1)*globalPoints_[2]
              +  0.25*(1.0 - eta1)*globalPoints_[3];

            J2D_[gp] = mag(dX_dEta1 ^ dX_dEta2);
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported face type for 2D Jacobian calculation: "
                << type_ << abort(FatalError);
        }
    }
}

// ************************************************************************* //
