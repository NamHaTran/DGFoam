/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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
    Foam::positivityPreservingLimiter

Description
    Positivity-preserving DG limiter following the 2D reference implementation.

\*---------------------------------------------------------------------------*/

#ifndef Foam_positivityPreservingLimiter_H
#define Foam_positivityPreservingLimiter_H

#include "dgLimiter.H"
#include "volFieldsFwd.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class positivityPreservingLimiter Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::positivityPreservingLimiter
 * \brief Zhang-Shu-style positivity limiter for DG conservative fields
 *
 * \details
 * This limiter enforces positivity by scaling higher-order modal content in
 * the conservative fields associated with density, momentum, and total energy.
 *
 * For each cell, the limiter:
 * - evaluates the minimum density over cell and face Gauss samples,
 * - computes \c theta1 for density positivity,
 * - evaluates the thermodynamic energy
 *   \f$\rho e = E - \tfrac{1}{2} |\rho U|^2/\rho\f$ on the modified state,
 * - computes \c theta2 for thermodynamic-energy positivity,
 * - scales higher-order modes of the limited fields by the required factor.
 *
 * The limiting pattern is:
 * \f[
 * (\rho,\rho U,E)
 * \;\rightarrow\;
 * \theta_1,\theta_2
 * \;\rightarrow\;
 * \text{scale higher-order modes}
 * \f]
 *
 * This implementation operates directly on every cell and therefore may be
 * used without a troubled-cell detector.
 */
class positivityPreservingLimiter
:
    public dgLimiter
{
    // Private data

        /** \brief Cached conservative state at one sampling point. */
        struct conservativeSample
        {
            scalar rho;
            vector rhoU;
            scalar E;
        };

        /** Name of the density field read from the limiter dictionary. */
        word densityFieldName_;

        /** Name of the momentum field read from the limiter dictionary. */
        word momentumFieldName_;

        /** Name of the total-energy field read from the limiter dictionary. */
        word energyFieldName_;

        /** Density field to be scaled by theta1*theta2. */
        dgField<scalar>* densityFieldPtr_;

        /** Momentum field used in the rhoe positivity check. */
        dgField<vector>* momentumFieldPtr_;

        /** Total-energy field used in the rhoe positivity check. */
        dgField<scalar>* energyFieldPtr_;

        /** Positivity tolerance read from dictionary keyword `tolerance`. */
        scalar epsilon_;

        /** Cell-local theta1 cached during correct(). */
        scalar theta1_;

        /** Cell-local theta2 cached during correct(). */
        scalar theta2_;

        /** Flag enabling runtime output of cell-wise theta coefficients. */
        bool writeTheta_;

        /** Optional runtime field storing the cell-wise theta1 coefficient. */
        mutable autoPtr<volScalarField> theta1FieldPtr_;

        /** Optional runtime field storing the cell-wise theta2 coefficient. */
        mutable autoPtr<volScalarField> theta2FieldPtr_;

        /** \brief Resolve density, momentum, and energy fields from the limiter setup. */
        void initializeConservativeFields();

        /** \brief Return the mean density stored in the zeroth modal coefficient. */
        scalar densityMean(const label cellID) const;

        /** \brief Return the theta1-modified density value at one sample point. */
        scalar modifiedDensityValue
        (
            const scalar rhoValue,
            const scalar meanRho,
            const scalar theta1
        ) const;

        /**
         * \brief Return \f$E - 0.5 |\rho U|^2/\rho\f$ with a guarded denominator.
         *
         * When the conserved scalar stores total energy, this corresponds to
         * \f$\rho e\f$. When it stores total enthalpy, it corresponds to
         * \f$\rho h\f$.
         */
        scalar calcRhoThermoEnergy
        (
            const scalar rho,
            const vector& rhoU,
            const scalar E
        ) const;

        /** \brief Clamp the Zhang-Shu style limiter coefficient to [0, 1]. */
        scalar thetaCoeff
        (
            const scalar meanValue,
            const scalar minValue,
            const scalar omega
        ) const;

        /**
         * \brief Compute theta1 and theta2 for one cell from a shared sample cache.
         *
         * The helper walks all cell/face Gauss samples once to gather
         * conservative states, then reuses that cache to evaluate the two
         * positivity coefficients.
         */
        void computeThetaCoeffs
        (
            const label cellID,
            scalar& theta1,
            scalar& theta2
        ) const;

        /** \brief Scale higher-order scalar modes and refresh the FOAM DoF view. */
        void scaleHighOrderModes
        (
            dgField<scalar>& field,
            const label cellID,
            const scalar theta
        ) const;

        /** \brief Scale higher-order vector modes and refresh the FOAM DoF view. */
        void scaleHighOrderModes
        (
            dgField<vector>& field,
            const label cellID,
            const scalar theta
        ) const;

        /** \brief Reset theta output fields to 1 on all cells. */
        void resetThetaFields() const;

        /** \brief Store one cell-local theta pair into the runtime fields. */
        void setThetaFields
        (
            const label cellID,
            const scalar theta1,
            const scalar theta2
        ) const;

protected:

        /** \brief Report limiter usage after one correct() pass. */
        virtual void postCorrect() override;

        /** \brief Limit one scalar field using the cell-local theta cache. */
        virtual void limitField
        (
            dgField<scalar>& field,
            const label cellID
        ) override;

        /** \brief Limit one vector field using theta2. */
        virtual void limitField
        (
            dgField<vector>& field,
            const label cellID
        ) override;

public:

    /** \brief Runtime type information. */
    TypeName("positivityPreserving");

    /**
     * \brief Construct the positivity-preserving limiter.
     *
     * Required dictionary entries:
     * \verbatim
     * density  rho;
     * momentum rhoU;
     * energy   E;
     * tolerance  1e-8;
     * writeTheta false;
     * \endverbatim
     *
     * Unlike detector-based limiters, this limiter may omit the base-class
     * \c detector entry because it evaluates positivity coefficients on every
     * cell directly inside \c correct().
     *
     * The optional base-class \c fields list is interpreted as additional
     * fields that should also be limited with \c theta2. The conservative
     * fields above are injected into the internal \c fields list before the
     * base-class constructor runs, so they are always limited and are looked
     * up only once by the base class.
     */
    positivityPreservingLimiter
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** \brief Re-read the limiter dictionary and refresh cached field pointers. */
    virtual void read(const dictionary& dict) override;

    /** \brief Compute theta1/theta2 on every cell and limit only cells needing correction. */
    virtual void correct() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
