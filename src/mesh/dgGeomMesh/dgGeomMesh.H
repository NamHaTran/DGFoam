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
    Foam::dgGeomMesh

Description
    Geometric mesh container for Discontinuous Galerkin discretization.

    Stores all DG geometric cells and faces constructed from an fvMesh
    and manages shared reference-cell and reference-face data.

SourceFiles
    dgGeomMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeomMesh_H
#define Foam_dgGeomMesh_H

#include "fvMesh.H"
#include "dgRefFace.H"
#include "dgGeomFace.H"
#include "dgGeomCell.H"
#include "dgRefCell.H"
#include "dgFacePosition.H"
#include "dgFaceType.H"
#include "dgCellType.H"

#include <memory>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgGeomMesh Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeomMesh
 * \brief DG geometric mesh wrapper around fvMesh
 *
 * \details
 * dgGeomMesh constructs and stores all DG-specific geometric
 * entities derived from an fvMesh.
 *
 * Responsibilities include:
 * - Creating and storing dgGeomFace and dgGeomCell objects
 * - Managing shared reference-face (dgRefFace) data
 * - Managing reference-cell (dgRefCell) data for different cell types
 * - Providing connectivity information required for DG surface and volume
 *   integration
 *
 * The class acts as the central geometric data container for
 * DG discretization.
 */
class dgGeomMesh
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Reference to the underlying fvMesh */
        const fvMesh& mesh_;

        /** Polynomial order for DG discretization */
        label pOrder_;

        /** Shared reference-face data */
        std::shared_ptr<dgRefFace> refFace_;

        /** Reference-cell data for tetrahedral cells */
        std::shared_ptr<dgRefCell> refCellTet_;

        /** Reference-cell data for hexahedral cells */
        std::shared_ptr<dgRefCell> refCellHex_;

        /** Reference-cell data for prismatic cells */
        std::shared_ptr<dgRefCell> refCellPrism_;

        /** Reference-cell data for pyramidal cells */
        std::shared_ptr<dgRefCell> refCellPyramid_;

        /** List of DG geometric faces */
        List<dgGeomFace*> faces_;

        /** List of DG geometric cells */
        List<dgGeomCell*> cells_;

        /** List of boundary face indices */
        labelList boundaryFaces_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct DG geometric mesh
         *
         * Builds all DG geometric faces and cells from the provided
         * fvMesh and initializes reference-cell and reference-face data
         * for the specified polynomial order.
         */
        dgGeomMesh(const fvMesh& mesh, const label pOrder);

    /** Destructor */
    ~dgGeomMesh();

    //====================================================================//
    //                         Access functions
    //====================================================================//

        /** Return list of DG geometric faces */
        inline const List<dgGeomFace*>& faces() const
        {
            return faces_;
        }

        /** Return list of DG geometric cells */
        inline const List<dgGeomCell*>& cells() const
        {
            return cells_;
        }

        /** Return DG polynomial order */
        inline label pOrder() const
        {
            return pOrder_;
        }

        /** Return number of DG faces */
        inline label nFaces() const
        {
            return faces_.size();
        }

        /** Return number of DG cells */
        inline label nCells() const
        {
            return cells_.size();
        }

        /** Return fvMesh face-owner list */
        inline const labelList& faceOwner() const
        {
            return mesh_.faceOwner();
        }

        /** Return number of internal faces */
        inline label nInternalFaces() const
        {
            return mesh_.nInternalFaces();
        }

    //====================================================================//
    //                        Boundary handling
    //====================================================================//

        /**
         * \brief Identify and store boundary faces
         */
        void getBoundaryFaces();

        /**
         * \brief Return list of boundary face indices
         */
        const Foam::labelList& boundaryFaces() const;

        /**
         * \brief Return patch ID for a given face
         */
        Foam::label getPatchID(const label faceID) const;

        /**
         * \brief Return local face index within a patch
         */
        Foam::label getLocalFaceID
        (
            const label faceID,
            const label patchID
        ) const;

    //====================================================================//
    //                        Mesh access
    //====================================================================//

        /** Return reference to the underlying fvMesh */
        inline const fvMesh& getFvMesh() const
        {
            return mesh_;
        }

    //====================================================================//
    //                    Parallel-processing support
    //====================================================================//

        /**
         * \brief Assign Gauss-point connectivity for processor-patch faces
         */
        void assignFaceConnectivityToProcPatches();
};

} // End namespace Foam

#endif

// ************************************************************************* //
