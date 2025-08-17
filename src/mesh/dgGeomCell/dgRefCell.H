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


#ifndef dgRefCell_H
#define dgRefCell_H

#include "label.H"
#include "scalar.H"
#include "List.H"
#include "vector.H"
#include "dgCellType.H"  // Enum defining HEX, PRISM, TET, etc.

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class dgRefCell Declaration
\*---------------------------------------------------------------------------*/

class dgRefCell
{
private:

    // Attributes

        // Polynomial order
        label pOrder_;

        // Cell type (e.g. HEX, PRISM, TET)
        dgCellType type_;

        // Number of Gauss points
        label nGauss_;

        // Degree of freedom (DOF) for the cell
        label nDof_;

        // Gauss points in reference cell: (eta1, eta2, eta3)
        List<vector> gaussP_;  // Each vector = (eta1, eta2, eta3)

        // Gauss weights
        List<scalar> wGauss_;

        // Basis functions and derivatives at Gauss points
        List<List<scalar>> basis_;
        List<List<scalar>> dBasis_dEta1_;
        List<List<scalar>> dBasis_dEta2_;
        List<List<scalar>> dBasis_dEta3_;

public:

    // Constructors

        //- Construct from order and cell type
        dgRefCell(const label pOrder, const dgCellType type);

    //- Destructor
    ~dgRefCell() = default;

    // Member Functions

        //- Return polynomial order
        inline label order() const { return pOrder_; }

        //- Return cell type
        inline dgCellType type() const { return type_; }

        //- Return number of Gauss points
        inline label nGauss() const { return nGauss_; }

        //- Return degree of freedom (DOF)
        inline label nDof() const { return nDof_; }

        //- Return Gauss points and weights
        inline const List<vector>& gaussPoints() const { return gaussP_; }
        inline const List<scalar>& weights() const     { return wGauss_; }

        //- Return basis and derivatives
        inline const List<List<scalar>>& basis() const         { return basis_; }
        inline const List<List<scalar>>& dBasis_dEta1() const  { return dBasis_dEta1_; }
        inline const List<List<scalar>>& dBasis_dEta2() const  { return dBasis_dEta2_; }
        inline const List<List<scalar>>& dBasis_dEta3() const  { return dBasis_dEta3_; }

        //- Generate Gauss points and weights for the given cell type
        void generateCellGaussPointsAndWeights();

        //- Evaluate basis and derivative functions at Gauss points
        void computeBasisAndDerivatives();
};
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
