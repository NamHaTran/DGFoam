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
    dgBasisFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgBasisField_H
#define Foam_dgBasisField_H

#include "GaussField.H"
#include "dgGeomMesh.H"
#include "Jacobian.H"

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

        /** \brief DG geometric mesh used to evaluate basis data. */
        const dgGeomMesh& mesh_;

        /** \brief Mesh cell index associated with this basis-field cache. */
        const label cellID_;

        /** \brief Local list of face labels belonging to the cell. */
        labelList faceLabels_;

        /** \brief Number of basis functions, equal to the local cell DoF count. */
        label nDof_;

        /** \brief Scalar basis-function values at cell and face Gauss points. */
        List<GaussField<scalar>> basis_;

        /** \brief Basis-function derivative values stored as vector Gauss fields. */
        List<GaussField<vector>> dBasis_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct basis-function storage for one DG cell.
         *
         * Initializes internal storage and precomputes the basis-function
         * fields and derivative fields associated with the selected cell.
         *
         * \param cellID Mesh cell index.
         * \param mesh DG geometric mesh that owns the cell.
         */
        dgBasisField(const label cellID, const dgGeomMesh& mesh);

    //====================================================================//
    //                        Basis computation
    //====================================================================//

        /**
         * \brief Compute scalar basis-function values at Gauss points.
         *
         * Populates one scalar GaussField per local basis function for the
         * cell interior and each face of the cell.
         */
        void computeBasisField();

        /**
         * \brief Compute derivative basis-function values at Gauss points.
         *
         * Populates one vector GaussField per local basis function using the
         * derivative data available for the associated cell and faces.
         */
        void computeDBasisField();

    //====================================================================//
    //                           Access
    //====================================================================//

        /**
         * \brief Return the scalar GaussField of one basis function.
         *
         * \param dof Local degree-of-freedom index.
         * \return Read-only basis-function values at Gauss points.
         */
        const GaussField<scalar>& getBasis(const label dof) const;

        /**
         * \brief Return the derivative GaussField of one basis function.
         *
         * \param dof Local degree-of-freedom index.
         * \return Read-only derivative values at Gauss points.
         */
        const GaussField<vector>& getDBasis(const label dof) const;

        /** \brief Return the owner cell index. */
        inline label cellID() const { return cellID_; }

        /** \brief Return the DG geometric mesh associated with this basis field. */
        inline const dgGeomMesh& dgMesh() const { return mesh_; }

        /** \brief Return the number of local degrees of freedom. */
        inline label nDof() const { return nDof_; }
};

} // End namespace Foam

// Include inline implementation
#include "dgBasisFieldI.H"

#endif

// ************************************************************************* //
