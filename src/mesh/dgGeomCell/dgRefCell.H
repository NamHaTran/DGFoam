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
    Foam::dgRefCell

Description
    Reference-cell definition for Discontinuous Galerkin discretization.

    Stores quadrature rules (Gauss points and weights) and basis-function
    evaluations in reference coordinates for a specific cell topology
    and polynomial order.

SourceFiles
    dgRefCell.C

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

/**
 * \class Foam::dgRefCell
 * \brief Reference-cell data for DG discretization
 *
 * \details
 * dgRefCell stores all reference-space data required by the
 * Discontinuous Galerkin method for a given cell type and
 * polynomial order.
 *
 * This includes:
 * - Gauss quadrature points and weights in reference coordinates
 * - Basis functions evaluated at Gauss points
 * - Derivatives of basis functions with respect to reference
 *   coordinates \f$\eta_1\f$, \f$\eta_2\f$, and \f$\eta_3\f$
 *
 * The data is generated once and shared among all physical cells
 * of the same type and order to avoid redundant computations.
 *
 * \ingroup grpDGMesh
 */
class dgRefCell
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Polynomial order of the DG approximation */
        label pOrder_;

        /** Reference cell type (HEX, PRISM, TET, etc.) */
        dgCellType type_;

        /** Number of Gauss quadrature points */
        label nGauss_;

        /** Number of degrees of freedom (DOF) */
        label nDof_;

        /** Gauss points in reference coordinates (eta1, eta2, eta3) */
        List<vector> gaussP_;

        /** Gauss quadrature weights */
        List<scalar> wGauss_;

        /** Basis functions evaluated at Gauss points */
        List<List<scalar>> basis_;

        /** Basis-function derivatives w.r.t. \f$\eta_1\f$ */
        List<List<scalar>> dBasis_dEta1_;

        /** Basis-function derivatives w.r.t. \f$\eta_2\f$ */
        List<List<scalar>> dBasis_dEta2_;

        /** Basis-function derivatives w.r.t. \f$\eta_3\f$ */
        List<List<scalar>> dBasis_dEta3_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct reference-cell data
         *
         * Initializes reference-space quadrature rules and basis
         * functions for the specified polynomial order and cell type.
         */
        dgRefCell(const label pOrder, const dgCellType type);

    /** Destructor */
    ~dgRefCell() = default;

    //====================================================================//
    //                         Access functions
    //====================================================================//

        /** Return polynomial order */
        inline label order() const { return pOrder_; }

        /** Return reference cell type */
        inline dgCellType type() const { return type_; }

        /** Return number of Gauss points */
        inline label nGauss() const { return nGauss_; }

        /** Return number of degrees of freedom */
        inline label nDof() const { return nDof_; }

        /** Return Gauss points in reference space */
        inline const List<vector>& gaussPoints() const { return gaussP_; }

        /** Return Gauss quadrature weights */
        inline const List<scalar>& weights() const     { return wGauss_; }

        /** Return basis functions evaluated at Gauss points */
        inline const List<List<scalar>>& basis() const
        {
            return basis_;
        }

        /** Return basis-function derivatives w.r.t. \f$\eta_1\f$ */
        inline const List<List<scalar>>& dBasis_dEta1() const
        {
            return dBasis_dEta1_;
        }

        /** Return basis-function derivatives w.r.t. \f$\eta_2\f$ */
        inline const List<List<scalar>>& dBasis_dEta2() const
        {
            return dBasis_dEta2_;
        }

        /** Return basis-function derivatives w.r.t. \f$\eta_3\f$ */
        inline const List<List<scalar>>& dBasis_dEta3() const
        {
            return dBasis_dEta3_;
        }

    //====================================================================//
    //                       Internal generation
    //====================================================================//

        /**
         * \brief Generate Gauss quadrature points and weights
         *
         * Generates reference-space quadrature rules based on
         * the cell type and polynomial order.
         */
        void generateCellGaussPointsAndWeights();

        /**
         * \brief Compute basis functions and derivatives
         *
         * Evaluates basis functions and their derivatives at
         * all Gauss quadrature points.
         */
        void computeBasisAndDerivatives();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
