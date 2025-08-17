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
    Geometric representation of a single cell in the mesh, tailored for use
    with the Discontinuous Galerkin (DG) method.

    The dgGeomCell class manages the geometric data of a finite volume mesh
    cell, including:
      - Its topological and geometric face list
      - Unique list of vertex coordinates
      - Access to polyMesh-based quantities such as volume and centroid

    During construction, the class initializes all dgGeomFace objects of the
    cell using a shared dgRefFace, and stores the geometry required for
    per-element numerical integration.

    The class also provides access to reference data (basis functions,
    Gauss points, weights) via a shared dgRefCell, which is used to evaluate
    solution fields within the cell.

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
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgGeomCell Declaration
\*---------------------------------------------------------------------------*/

class dgGeomCell
{
private:

    //- Cell ID in the mesh
    label cellID_;

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Cell type
    dgCellType type_;

    //- Shared pointer to reference cell
    std::shared_ptr<dgRefCell> refCell_;

    //- Shared pointer to reference face (used to create dgGeomFace)
    std::shared_ptr<dgRefFace> refFace_;

    //- List of face ids of the cell
    labelList faceLabels_;

    //- List of unique vertices of the cell
    pointField cellPoints_;

    //- List of point labels in the cell (for ordered access)
    List<label> cellPointLabels_;

public:

    // Constructors

        //- Construct from cell ID, mesh and reference cell/face
        dgGeomCell
        (
            label cellID,
            const fvMesh& mesh,
            const std::shared_ptr<dgRefCell>& refCell
        );

    //- Destructor
    ~dgGeomCell();

    // Member Functions

        //- Return the list of face ids of the cell
        const labelList& faces() const;

        //- Return the number of faces in the cell
        label nFaces() const;

        //- Return list of unique cell points
        const pointField& points() const;

        //- Return the number of unique points in the cell
        label nPoints() const;

        //- Return cell centroid (from mesh)
        point centre() const;

        //- Return cell volume (from mesh)
        scalar volume() const;

        //- Update face information based on the cell's geometry and mesh
        void updateFaceInfo
        (
            List<dgGeomFace*>& faces
        );

    // Access to reference cell data

        //- Return Gauss points in reference cell
        inline const List<vector>& gaussPoints() const
        {
            return refCell_->gaussPoints();
        }

        //- Return number of Gauss points in reference cell
        inline scalar nGauss() const
        {
            return refCell_->nGauss();
        }

        //- Return weights of Gauss points
        inline const List<scalar>& weights() const
        {
            return refCell_->weights();
        }

        //- Return basis functions at Gauss points
        inline const List<List<scalar>>& basis() const
        {
            return refCell_->basis();
        }

        //- Return number of Dof in reference cell
        inline scalar nDof() const
        {
            return refCell_->nDof();
        }

        //- Return derivative of basis functions w.r.t η₁
        inline const List<List<scalar>>& dBasis_dEta1() const
        {
            return refCell_->dBasis_dEta1();
        }

        //- Return derivative of basis functions w.r.t η₂
        inline const List<List<scalar>>& dBasis_dEta2() const
        {
            return refCell_->dBasis_dEta2();
        }

        //- Return derivative of basis functions w.r.t η₃
        inline const List<List<scalar>>& dBasis_dEta3() const
        {
            return refCell_->dBasis_dEta3();
        }

        //- Print debug information about the cell
        void printDebugInfo() const;
};

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
