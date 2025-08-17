/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Class
    Foam::dgGeomFace

Description
    Geometric face of DG.

See also

SourceFiles
    dgGeomFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeomFace_H
#define Foam_dgGeomFace_H

#include "face.H"
#include "fvMesh.H"
#include "pointField.H"
#include "scalarField.H"
#include "dgRefFace.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "dgFacePosition.H"
#include "dgFaceType.H"
#include "dgCellType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class Foam::dgRefFace;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgGeomFace Declaration
\*---------------------------------------------------------------------------*/

class dgGeomFace
{
private:

    // Attributes

        //- face ID in the mesh
        label faceID_;

        //- Fave type
        dgFaceType type_;

        //- Reference to the face (list of point indices)
        const fvMesh& mesh_; 

        //- Shared reference face data (precomputed basis, Gauss points, weights)
        std::shared_ptr<dgRefFace> refFace_;  
        
        //- Points list
        List<vector> globalPoints_;
        
        //- Face position in owner
        dgFacePosition ownerPos_;

        //- Face position in neighbor
        dgFacePosition neighborPos_;

        //- Type of owner cell
        dgCellType ownerCellType_;

        //- Type of neighbor cell
        dgCellType neighborCellType_;

        //- Global point IDs sorted in counter-clockwise order
        // List<label> sortedPointLabels_;

        //- Flattened 2D coordinates of the points
        List<vector2D> flattenedPoints_;

        //- Gauss point connectivity of owner and neighbor
        List<label> connectivity_;

public:

    // Constructors

        //- Default constructor
        dgGeomFace();

        //- Constructor from face and mesh
        dgGeomFace
        (
            label faceID,
            const fvMesh& mesh,
            const std::shared_ptr<dgRefFace>& refFace
        );

        //- Copy constructor
        dgGeomFace(const dgGeomFace& other);

    //- Destructor
    ~dgGeomFace() = default;

    // Member Functions
        //- Return number of points in the face
        label size() const;

        //- Return reference to point labels (from face)
        const face& baseFace() const;

        //- Return coordinates of the i-th point on the face
        const point& getPoint(label i) const;

        //- Return face center (calls face::centre internally)
        point centre() const;

        //- Return normal vector of the face (flat face assumed)
        vector normal() const;

        //- Return area normal vector of the face
        vector areaNormal() const;

        //- Return face area
        scalar area() const;

        //- Return face ID
        inline const label& id() const
        {
            return faceID_;
        }

        //- Flatten face points to 2D using face normal
        void flattenFace();

        //- Sort flattened face points in counter-clockwise order
        // void sortPointsCCW();

        //- Perform flattening, sorting and store data
        void processFlatAndSortedPoints();

        //- Set face position in owner
        inline void setOwnerPos(const dgFacePosition pos)
        {
            ownerPos_ = pos;
        }

        //- Set face position in neighbor
        inline void setNeighborPos(const dgFacePosition pos)
        {
            neighborPos_ = pos;
        }

        //- Set owner cell type
        inline void setOwnerCellType(const dgCellType type)
        {
            ownerCellType_ = type;
        }

        //- Set neighbor cell type
        inline void setNeighborCellType(const dgCellType type)
        {
            neighborCellType_ = type;
        }

        // Compute connectivity between Gauss points from owner and neighbor
        void findGaussConnectivity();

        // Mapping from reference to physical space for quad or tri face
        void mappingFromRefToReal
        (
            const dgFaceType type,
            const List<vector>& gaussPoints,
            const List<vector>& faceVertices,
            List<vector2D>& physicGaussP
        );
    
    // Access to reference face data

        //- Return Gauss points in reference face base on position on the owner
        inline const List<vector>& gaussPointsOwner() const
        {
            return refFace_->gaussPoints(ownerPos_);
        }

        //- Return Gauss points in reference face base on position on the neighbor
        inline const List<vector>& gaussPointsNeighbor() const
        {
            return refFace_->gaussPoints(neighborPos_);
        }

        //- Return weights of Gauss points
        inline const List<scalar>& weights() const
        {
            return refFace_->weights();
        }

        //- Return sorted point labels in counter-clockwise order
        /*
        inline const List<label>& sortedLabels() const
        {
            return sortedPointLabels_;
        }
        */

        //- Return flattened 2D coordinates of the points
        inline const List<vector2D>& flatCoords() const
        {
            return flattenedPoints_;
        }

        //- Return face owner list from fvMesh
        inline const labelList& faceOwner() const { return mesh_.faceOwner(); }

        //- Return number of internal faces from fvMesh
        inline label nInternalFaces() const { return mesh_.nInternalFaces(); }

        //- Return face connectivity
        inline const List<label>& connectivity() const
        {
            return connectivity_;
        }

    //- Operator overload
    
        // Copy assignment
        dgGeomFace& operator=(const dgGeomFace& other);

    //- Print debug information about the face
        void printDebugInfo() const;
};
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
