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

Class
    Foam::dgGeomFace

Description
    Geometric representation of a DG face.

    This class stores geometric and topological information associated
    with a mesh face for Discontinuous Galerkin (DG) discretization,
    including face orientation, reference-face mapping, Gauss-point
    connectivity, and basis-function data on both owner and neighbour
    sides.

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
#include "dgGeomFaceTools.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "dgFacePosition.H"
#include "dgFaceType.H"
#include "dgCellType.H"
#include "Jacobian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgGeomFace Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeomFace
 * \brief Geometric face data for DG discretization
 *
 * \details
 * dgGeomFace represents a single mesh face in the DG framework.
 * It provides:
 * - Access to geometric quantities (area, normal, centre)
 * - Face orientation relative to owner and neighbour cells
 * - Mapping between reference and physical face coordinates
 * - Gauss-point connectivity and basis-function data for surface integrals
 *
 * The face uses shared reference-face data (dgRefFace) to evaluate
 * basis functions and quadrature rules consistently across all faces
 * of the same type.
 */
class dgGeomFace
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Face index in the fvMesh */
        label faceID_;

        /** DG face type (TRI or QUAD) */
        dgFaceType type_;

        /** Reference to the finite-volume mesh */
        const fvMesh& mesh_;

        /** Shared reference-face data (basis, Gauss points, weights) */
        std::shared_ptr<dgRefFace> refFace_;

        /** Number of Gauss points on face */
        label nGauss_;

        /** Global face vertex coordinates */
        List<vector> globalPoints_;

        /** Owner cell ID */
        label ownerId_;

        /** Neighbor cell ID (-1 for boundary faces) */
        label neighborId_;

        /** Face position in the owner cell reference element */
        dgFacePosition ownerPos_;

        /** Face position in the neighbour cell reference element */
        dgFacePosition neighborPos_;

        /** Owner cell type */
        dgCellType ownerCellType_;

        /** Neighbour cell type */
        dgCellType neighborCellType_;

        /** Flattened 2D coordinates of face vertices */
        List<vector2D> flattenedPoints_;

        /** Physical coordinates of the face Gauss points */
        List<vector> physGaussPoints_;

        /** Gauss points on the face, respect to owner cell */
        List<vector> gaussPointsOwner_;

        /** Gauss points on the face, respect to neighbour cell */
        List<vector> gaussPointsNeighbor_;

        /** Basis-function data on the owner side */
        Foam::basisData ownerBasisData_;

        /** Basis-function data on the neighbour side */
        Foam::basisData neighborBasisData_;

        /** Jacobian 2D of faces */
        List<scalar> J2D_;

        //====================== Boundary information ======================//

            /** Patch index owning this face (-1 for internal faces) */
            label patchID_;

            /** True if this face lies on a boundary */
            bool isBoundary_;

            /** True if this face lies on a processor patch */
            bool isProcessorPatch_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** Default constructor */
        dgGeomFace();

        /**
         * \brief Construct from face ID, mesh, and reference-face data
         */
        dgGeomFace
        (
            label faceID,
            const fvMesh& mesh,
            const std::shared_ptr<dgRefFace>& refFace
        );

        /** Copy constructor */
        dgGeomFace(const dgGeomFace& other);

    /** Destructor */
    ~dgGeomFace() = default;

    //====================================================================//
    //                         Geometric access
    //====================================================================//

        /** Return number of vertices of the face */
        label size() const;

        /** Return nunber of Gauss points on the face */
        inline label nGauss() const { return nGauss_; }

        /** Return underlying OpenFOAM face */
        const face& baseFace() const;

        /** Return coordinate of the i-th face vertex */
        const point& getPoint(label i) const;

        /** Return face centroid */
        point centre() const;

        /** Return face unit normal (planar assumption) */
        vector normal() const;

        /** Return area-normal vector */
        vector areaNormal() const;

        /** Return face area */
        scalar area() const;

        /** Return face index */
        inline const label& id() const { return faceID_; }

        /** Return Jacobian 2D at Gauss points */
        inline const List<scalar>& J2D() const
        {
            return J2D_;
        }

    //====================================================================//
    //                       Face processing utilities
    //====================================================================//

        /** \brief Compute Gauss-point connectivity between owner and neighbour cells
         * 
         *  \param globalPoints Global coordinates of the face vertices (in the order of owner cell)
         *  \param flattenedPoints Flattened 2D coordinates of the face vertices (in the order of owner cell)
         */
        void flattenFace
        (
            const List<vector>& globalPoints,   // Input 3D coordinates
            List<vector2D>& flattenedPoints     // Output 2D coordinates
        ) const;

        /**
         * \brief Perform face flattening and ordering
         *
         * Computes flattened coordinates and prepares data required
         * for basis-function evaluation on the face.
         */
        void processFlatAndSortedPoints();

        /** Set face owner cell ID */
        inline void setOwnerId(const label id) { ownerId_ = id; }

        /** Set face neighbor cell ID */
        inline void setNeighborId(const label id) { neighborId_ = id; }

        /** Set face position in owner cell */
        inline void setOwnerPos(const dgFacePosition pos) { ownerPos_ = pos; }

        /** Set face position in neighbour cell */
        inline void setNeighborPos(const dgFacePosition pos) { neighborPos_ = pos; }

        /** Set owner cell type */
        inline void setOwnerCellType(const dgCellType type) { ownerCellType_ = type; }

        /** Set neighbour cell type */
        inline void setNeighborCellType(const dgCellType type) { neighborCellType_ = type; }

        /** Compute Gauss-point connectivity between owner and neighbour */
        void findGaussConnectivity();

        /** Compute basis functions and derivatives on the face */
        void computeBasisAndDerivatives();

        /**
         * \brief Compute Lame parameters at face Gauss points from both owner and neighbour cells
         * This method MUST be called after setting ownerPos_, neighborPos_, ownerCellType_ and neighborCellType_.
         */
        void computeLameParameters();

    //====================================================================//
    //                          Boundary queries
    //====================================================================//

        /** Return true if this face is a boundary face */
        inline bool isBoundary() const { return isBoundary_; }

        /** Return patch index of the face */
        inline label patchID() const { return patchID_; }

        /** Return true if this face lies on a processor patch */
        inline bool isProcessorPatch() const { return isProcessorPatch_; }

        /** Check if the given cell is the owner of this face */
        inline bool isOwner(const label cellID) const
        {
            return (cellID == mesh_.faceOwner()[faceID_]);
        }

    //====================================================================//
    //                    Reference-face data access
    //====================================================================//
        /** Return physical Gauss points */
        inline const List<vector>& physicGaussPoints() const
        {
            return physGaussPoints_;
        }

        /** Return Gauss points on owner-side reference face */
        inline const List<vector>& gaussPointsOwner() const
        {
            return gaussPointsOwner_;
        }

        /** Return Gauss points on neighbour-side reference face */
        inline const List<vector>& gaussPointsNeighbor() const
        {
            return gaussPointsNeighbor_;
        }

        /** Return Gauss quadrature weights */
        inline const List<scalar>& weights() const
        {
            return refFace_->weights();
        }

        /** Return flattened 2D face coordinates */
        inline const List<vector2D>& flatCoords() const
        {
            return flattenedPoints_;
        }

        /** Return fvMesh face-owner list */
        inline const labelList& faceOwner() const { return mesh_.faceOwner(); }

        /** Return number of internal faces */
        inline label nInternalFaces() const { return mesh_.nInternalFaces(); }

        /** Return owner-side basis functions */
        inline const List<List<scalar>>& ownerBasis() const
        {
            return ownerBasisData_.basis;
        }

        /** Return neighbour-side basis functions */
        inline const List<List<scalar>>& neighborBasis() const
        {
            return neighborBasisData_.basis;
        }

        /** Return owner-side basis derivatives */
        inline const List<List<scalar>>& owner_dBasis_dEta1() const
        {
            return ownerBasisData_.dBasis_dEta1;
        }
        inline const List<List<scalar>>& owner_dBasis_dEta2() const
        {
            return ownerBasisData_.dBasis_dEta2;
        }
        inline const List<List<scalar>>& owner_dBasis_dEta3() const
        {
            return ownerBasisData_.dBasis_dEta3;
        }

        /** Return neighbour-side basis derivatives */
        inline const List<List<scalar>>& neighbor_dBasis_dEta1() const
        {
            return neighborBasisData_.dBasis_dEta1;
        }
        inline const List<List<scalar>>& neighbor_dBasis_dEta2() const
        {
            return neighborBasisData_.dBasis_dEta2;
        }
        inline const List<List<scalar>>& neighbor_dBasis_dEta3() const
        {
            return neighborBasisData_.dBasis_dEta3;
        }

    //====================================================================//
    //                        Parallel processing
    //====================================================================//

        /** Mark this face as a processor-patch face */
        inline void setProcessorPatch(const bool isProcPatch)
        {
            isProcessorPatch_ = isProcPatch;
        }

        /** Set Gauss points */
        inline void setPhysicGaussPoints(const List<vector>& points)
        {
            physGaussPoints_ = points;
        }

        /** Set Lame parameters */
        inline void setLameParams(const List<scalar>& params)
        {
            J2D_ = params;
        }

        /** Calculate Gauss points from owner side, bases on physic Gauss */
        void calcGaussPointsOnOwnerSide();

    //====================================================================//
    //                          Operators & debug
    //====================================================================//

        /** Copy assignment operator */
        dgGeomFace& operator=(const dgGeomFace& other);

        /** Print debug information */
        void printDebugInfo() const;

        /** Print owner and neighbour cell information */
        void printCellsInfo() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
