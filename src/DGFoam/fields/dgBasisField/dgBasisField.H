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
    Foam::dgBasisField

Description
    Basis-function Gauss-field container for DG discretization.

    Stores precomputed GaussField representations of all basis functions
    and their derivatives for a single DG cell.

SourceFiles
    dgBasisField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgBasisField_H
#define Foam_dgBasisField_H

#include "GaussField.H"
#include "dgGeomMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dgBasisField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgBasisField
 * \brief Precomputed basis-function fields on Gauss points
 *
 * \details
 * dgBasisField stores GaussField representations of all polynomial
 * basis functions and their spatial derivatives associated with a
 * single DG cell.
 *
 * The basis functions are evaluated at Gauss points inside the cell
 * and on its faces using geometric information provided by dgGeomMesh.
 *
 * This class is typically used during DG residual assembly to evaluate
 * volume and surface integrals efficiently without recomputing basis
 * functions at runtime.
 */
class dgBasisField
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Reference to the DG geometric mesh */
        const dgGeomMesh& mesh_;

        /** Cell index associated with this basis field */
        const label cellID_;

        /** Local face labels of the cell */
        labelList faceLabels_;

        /** Number of degrees of freedom in this cell */
        label nDof_;

        /** Scalar basis-function Gauss fields (one per DoF) */
        List<GaussField<scalar>> basis_;

        /** Derivative basis-function Gauss fields (one per DoF) */
        List<GaussField<vector>> dBasis_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct basis fields for a given cell
         *
         * Initializes internal storage based on the cell geometry
         * and polynomial order defined in the DG mesh.
         */
        dgBasisField(const label cellID, const dgGeomMesh& mesh);

    //====================================================================//
    //                        Basis computation
    //====================================================================//

        /**
         * \brief Compute scalar basis-function Gauss fields
         *
         * Computes GaussField representations of all scalar
         * basis functions for the associated cell.
         */
        void computeBasisField();

        /**
         * \brief Compute derivative basis-function Gauss fields
         *
         * Computes GaussField representations of spatial derivatives
         * of all basis functions for the associated cell.
         */
        void computeDBasisField();

    //====================================================================//
    //                           Access
    //====================================================================//

        /**
         * \brief Return GaussField of the specified basis function
         *
         * \param dof Degree-of-freedom index
         */
        const GaussField<scalar>& getBasis(const label dof) const;

        /**
         * \brief Return GaussField of the specified derivative basis function
         *
         * \param dof Degree-of-freedom index
         */
        const GaussField<vector>& getDBasis(const label dof) const;

        /** Return number of degrees of freedom */
        inline label nDof() const { return nDof_; }
};

} // End namespace Foam

// Include inline implementation
#include "dgBasisFieldI.H"

#endif

// ************************************************************************* //
