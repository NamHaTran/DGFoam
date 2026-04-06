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
    Foam::meanValueLimiter

Description
    Mean-value limiter for conservative DG states.

\*---------------------------------------------------------------------------*/

#ifndef Foam_meanValueLimiter_H
#define Foam_meanValueLimiter_H

#include "dgLimiter.H"
#include "dgThermoConservative.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meanValueLimiter Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::meanValueLimiter
 * \brief Mean-state clamp for DG conservative variables
 *
 * \details
 * This limiter modifies only the zeroth modal coefficient of the conservative
 * state in each cell and resets selected higher-order modes when a clamp is
 * triggered.
 *
 * For each cell, the limiter applies the following sequence:
 * - clamp \f$\rho_0\f$ into \f$[\rho_{\min}, \rho_{\max}]\f$,
 * - compute the mean velocity from the stored unclamped density and clamp
 *   \f$|\mathbf{U}_0|\f$ by \c magUMax,
 * - reconstruct the mean temperature from \f$(\rho_0,\rho U_0,E_0)\f$,
 * - clamp \f$T_0\f$ into \f$[T_{\min}, T_{\max}]\f$ and rebuild \f$E_0\f$
 *   through \c dgThermoConservative.
 *
 * Whenever a variable is clamped, the corresponding higher-order modes are
 * set to zero so the updated mean state remains consistent with the retained
 * modal representation.
 *
 * The limiter acts on every cell directly and does not use a troubled-cell
 * detector.
 */
class meanValueLimiter
:
    public dgLimiter
{
    // Private data

        /** Name of the density field read from the limiter dictionary. */
        word densityFieldName_;

        /** Name of the momentum field read from the limiter dictionary. */
        word momentumFieldName_;

        /** Name of the total-energy field read from the limiter dictionary. */
        word energyFieldName_;

        /** Density field whose cell mean may be clamped. */
        dgField<scalar>* densityFieldPtr_;

        /** Momentum field whose cell mean may be clamped. */
        dgField<vector>* momentumFieldPtr_;

        /** Total-energy field whose cell mean may be clamped. */
        dgField<scalar>* energyFieldPtr_;

        /** Conservative thermodynamics used for T/E conversions. */
        const dgThermoConservative* thermoPtr_;

        /** Guard used when dividing by a stored unclamped density close to zero. */
        scalar epsilon_;

        /** Minimum allowed density mean value. */
        scalar rhoMin_;

        /** Maximum allowed density mean value. */
        scalar rhoMax_;

        /** Maximum allowed mean velocity magnitude. */
        scalar magUMax_;

        /** Minimum allowed mean temperature. */
        scalar Tmin_;

        /** Maximum allowed mean temperature. */
        scalar Tmax_;

        /** \brief Resolve density, momentum, and energy fields from the limiter setup. */
        void initializeConservativeFields();

        /** \brief Resolve the conservative thermodynamic model from the mesh registry. */
        void initializeThermo();

        /** \brief Return a sign-preserving safe density for division. */
        scalar safeDensityForDivision(const scalar rho) const;

        /** \brief Return the mean-state temperature reconstructed from conservative variables. */
        scalar calcMeanTemperature
        (
            const scalar rho,
            const vector& rhoU,
            const scalar E
        ) const;

        /** \brief Return the conservative total energy reconstructed from mean primitive variables. */
        scalar calcMeanTotalEnergy
        (
            const scalar rho,
            const vector& U,
            const scalar T
        ) const;

        /** \brief Set all higher-order scalar modes in one cell to zero. */
        void zeroHighOrderModes
        (
            dgField<scalar>& field,
            const label cellID
        ) const;

        /** \brief Set all higher-order vector modes in one cell to zero. */
        void zeroHighOrderModes
        (
            dgField<vector>& field,
            const label cellID
        ) const;

        /** \brief Clamp the conservative mean state of one cell and return true when modified. */
        bool limitCellMeanValue(const label cellID);

protected:

        /** \brief Report limiter usage after one correct() pass. */
        virtual void postCorrect() override;

        /** \brief Mean-value limiting is handled centrally in correct(). */
        virtual void limitField
        (
            dgField<scalar>& field,
            const label cellID
        ) override;

        /** \brief Mean-value limiting is handled centrally in correct(). */
        virtual void limitField
        (
            dgField<vector>& field,
            const label cellID
        ) override;

public:

    /** \brief Runtime type information. */
    TypeName("meanValue");

    /**
     * \brief Construct the mean-value limiter.
     *
     * Required dictionary entries:
     * \verbatim
     * density  rho;
     * momentum rhoU;
     * energy   E;
     * rhoMin   1e-8;
     * rhoMax   1e10;
     * magUMax  1e10;
     * Tmin     1e-8;
     * Tmax     1e10;
     * \endverbatim
     *
     * The limiter operates on every cell directly inside \c correct() and
     * therefore does not require a trouble-cell detector.
     */
    meanValueLimiter
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** \brief Re-read the limiter dictionary and refresh cached field pointers. */
    virtual void read(const dictionary& dict) override;

    /** \brief Clamp rho/U/T mean values on every cell. */
    virtual void correct() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
