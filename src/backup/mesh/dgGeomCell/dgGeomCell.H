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
    Foam::dgGeomCell

Description
    Geometric representation of a single DG cell.

    Stores topological and geometric information required for
    Discontinuous Galerkin (DG) discretization, including connectivity,
    vertex coordinates, and access to reference-cell data for numerical
    integration.

SourceFiles
    dgGeomCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeomCell_H
#define Foam_dgGeomCell_H

#include "fvMesh.H"
#include "pointField.H"
#include "scalarField.H"
#include "dgRefCell.H"
#include "dgRefFace.H"
#include "dgGeomFace.H"
#include "volFields.H"
#include "Jacobian.H"
#include "SquareMatrix.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgGeomCell Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeomCell
 * \brief Geometric container for a DG cell
 *
 * \details
 * The dgGeomCell class represents the geometry of a single mesh cell
 * tailored for Discontinuous Galerkin discretization.
 *
 * It provides:
 * - Access to mesh-based quantities (volume, centroid)
 * - Cell topology (faces, neighbors, vertices)
 * - Access to reference-cell data (Gauss points, weights, basis functions)
 *
 * Geometry-specific face objects (dgGeomFace) are initialized using a
 * shared reference face (dgRefFace). Reference-cell data is shared via
 * dgRefCell to avoid duplication and ensure consistency.
 *
 * \ingroup grpDGMesh
 */
class dgGeomCell
{
private:

    //====================================================================//
    //                        Data members
    //====================================================================//

    /** Cell index in the mesh */
    label cellID_;

    /** Reference to the finite-volume mesh */
    const fvMesh& mesh_;

    /** DG cell type */
    dgCellType type_;

    /** Shared reference-cell data (basis, Gauss points, weights) */
    std::shared_ptr<dgRefCell> refCell_;

    /** Shared reference-face data used to construct dgGeomFace */
    std::shared_ptr<dgRefFace> refFace_;

    /** List of face indices belonging to this cell */
    labelList faceLabels_;

    /** List of neighbouring cell indices */
    labelList neighborCellLabels_;

    /** Unique list of vertex coordinates of the cell */
    pointField cellPoints_;

    /** Ordered list of vertex indices for consistent access */
    List<label> cellPointLabels_;

    /** List of internal Gauss Jacobian */
    List<scalar> J3D_;

    /** Mass matrix for the cell */
    Foam::SquareMatrix<scalar> massMatrix_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct from mesh cell ID and reference data
         *
         * Initializes geometric connectivity and associates the cell
         * with shared reference-cell and reference-face data.
         */
        dgGeomCell
        (
            label cellID,
            const fvMesh& mesh,
            const std::shared_ptr<dgRefCell>& refCell
        );

    /** Destructor */
    ~dgGeomCell();

    //====================================================================//
    //                        Topological access
    //====================================================================//

        /** Return list of face indices of the cell */
        inline const labelList& faces() const
        {
            return faceLabels_;
        }

        /** Return list of neighbouring cell indices */
        inline const labelList& neighborCells() const
        {
            return neighborCellLabels_;
        }

        /** Return number of faces of the cell */
        inline label nFaces() const
        {
            return faceLabels_.size();
        }

        /** Return unique vertex coordinates of the cell */
        inline const pointField& points() const
        {
            return cellPoints_;
        }

        /** Return number of unique vertices */
        inline label nPoints() const
        {
            return cellPoints_.size();
        }

    //====================================================================//
    //                    Mesh-based geometric data
    //====================================================================//

        /** Return cell centroid (from fvMesh) */
        point centre() const;

        /** Return cell volume (from fvMesh) */
        scalar volume() const;

        /**
         * \brief Update face geometry information
         *
         * Initializes or updates dgGeomFace objects associated with
         * this cell based on mesh connectivity.
         */
        void updateFaceInfo
        (
            List<dgGeomFace*>& faces
        );

        /** Return the mass matrix for the cell */
        inline const Foam::SquareMatrix<scalar>& massMatrix() const
        {
            return massMatrix_;
        }

    //====================================================================//
    //                 Reference-cell numerical data
    //====================================================================//

        /** Return Gauss points in the reference cell */
        inline const List<vector>& gaussPoints() const
        {
            return refCell_->gaussPoints();
        }

        /** Return number of Gauss points */
        inline scalar nGauss() const
        {
            return refCell_->nGauss();
        }

        /** Return Gauss weights */
        inline const List<scalar>& weights() const
        {
            return refCell_->weights();
        }

        /** Return basis functions evaluated at Gauss points */
        inline const List<List<scalar>>& basis() const
        {
            return refCell_->basis();
        }

        /** Return number of degrees of freedom */
        inline scalar nDof() const
        {
            return refCell_->nDof();
        }

        /** Return derivative of basis functions w.r.t. \f$\eta_1\f$ */
        inline const List<List<scalar>>& dBasis_dEta1() const
        {
            return refCell_->dBasis_dEta1();
        }

        /** Return derivative of basis functions w.r.t. \f$\eta_2\f$ */
        inline const List<List<scalar>>& dBasis_dEta2() const
        {
            return refCell_->dBasis_dEta2();
        }

        /** Return derivative of basis functions w.r.t. \f$\eta_3\f$ */
        inline const List<List<scalar>>& dBasis_dEta3() const
        {
            return refCell_->dBasis_dEta3();
        }

        /** Return Jacobian 3D at Gauss points */
        inline const List<scalar>& J3D() const
        {
            return J3D_;
        }

    //====================================================================//
    //                           Utilities
    //====================================================================//

        /** Print debug information about the cell */
        void printDebugInfo() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
