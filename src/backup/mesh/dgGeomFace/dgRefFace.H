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
    Foam::dgRefFace

Description
    Reference-face definition for DG discretization.

    Provides Gauss quadrature points, weights, and basis-function data
    on reference faces of 3D elements following standard DG reference
    element conventions.

SourceFiles
    dgRefFace.C

\*---------------------------------------------------------------------------*/

#ifndef dgRefFace_H
#define dgRefFace_H

#include "label.H"
#include "scalar.H"
#include "List.H"
#include "vector.H"
#include "dgFacePosition.H"
#include "dgCellType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Struct basisData Declaration
\*---------------------------------------------------------------------------*/

/**
 * \struct basisData
 * \brief Basis-function data on a reference face
 *
 * \details
 * Stores basis-function values and their derivatives evaluated
 * at Gauss quadrature points on a reference face.
 *
 * The derivatives are expressed with respect to the reference
 * coordinates \f$\eta_1\f$, \f$\eta_2\f$, and \f$\eta_3\f$.
 */
struct basisData
{
    /** Basis-function values at Gauss points */
    List<List<scalar>> basis;

    /** Derivatives w.r.t. \f$\eta_1\f$ */
    List<List<scalar>> dBasis_dEta1;

    /** Derivatives w.r.t. \f$\eta_2\f$ */
    List<List<scalar>> dBasis_dEta2;

    /** Derivatives w.r.t. \f$\eta_3\f$ */
    List<List<scalar>> dBasis_dEta3;
};


/*---------------------------------------------------------------------------*\
                           Class dgRefFace Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgRefFace
 * \brief Reference-face data for DG surface integration
 *
 * \details
 * dgRefFace represents a reference face element used in
 * Discontinuous Galerkin discretization of 3D meshes.
 *
 * The class provides:
 * - Gauss quadrature points and weights on reference faces
 * - Face-position-dependent mappings in the reference cube
 * - Basis functions and their derivatives evaluated on faces
 *
 * The reference faces follow the vertex ordering and orientation
 * conventions described in:
 *
 * Sherwin, S.J., Karniadakis, G.E.,
 * "Spectral/hp Element Methods for Computational Fluid Dynamics",
 * 2nd Edition, Oxford University Press, 2005.
 *
 * Data stored in this class is shared among all geometric faces
 * of the same polynomial order to avoid redundant computation.
 */
class dgRefFace
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Polynomial order of the DG approximation */
        label pOrder_;

        /** Number of Gauss quadrature points */
        label nGauss_;

        /** Gauss quadrature weights (shared for all faces) */
        List<scalar> wGauss_;

        /** Gauss points on face ABCD in reference cube */
        List<vector> gaussP_ABCD_;

        /** Gauss points on face EFGH in reference cube */
        List<vector> gaussP_EFGH_;

        /** Gauss points on face ABEF in reference cube */
        List<vector> gaussP_ABEF_;

        /** Gauss points on face CDGH in reference cube */
        List<vector> gaussP_CDGH_;

        /** Gauss points on face BDFH in reference cube */
        List<vector> gaussP_BDFH_;

        /** Gauss points on face ACEG in reference cube */
        List<vector> gaussP_ACEG_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct reference-face data
         *
         * Initializes Gauss quadrature rules for the specified
         * polynomial order.
         */
        dgRefFace(const label pOrder);

    /** Destructor */
    ~dgRefFace() = default;

    //====================================================================//
    //                         Access functions
    //====================================================================//

        /** Return polynomial order */
        inline label order() const { return pOrder_; }

        /** Return number of Gauss points */
        inline label nGauss() const { return nGauss_; }

        /** Return Gauss quadrature weights */
        inline const List<scalar>& weights() const { return wGauss_; }

        /**
         * \brief Return Gauss points for a given face position
         *
         * The returned Gauss points are expressed in the
         * reference-cube coordinate system.
         */
        const List<vector>& gaussPoints(const dgFacePosition facePos) const;

    //====================================================================//
    //                      Internal generation
    //====================================================================//

        /**
         * \brief Generate Gauss points and weights for all face positions
         *
         * Generates quadrature rules on reference faces based on
         * the polynomial order.
         */
        void generateFaceGaussPointsAndWeights();

        /**
         * \brief Compute basis functions and derivatives on a reference face
         *
         * Computes basis-function values and their derivatives
         * for the specified cell type and face position.
         */
        Foam::basisData computeBasisAndDerivatives
        (
            const dgCellType cellType,
            const dgFacePosition pos
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
