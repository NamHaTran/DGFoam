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

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgArtificialDiffusionNonSmooth_H
#define Foam_dgArtificialDiffusionNonSmooth_H

#include "IOdictionary.H"
#include "dgField.H"
#include "dgCellType.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dgArtificialDiffusionNonSmooth Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgArtificialDiffusionNonSmooth
 * \brief NonSmooth artificial diffusion coefficient.
 *
 * Reads \c artificialDiffusion from \c constant/dgOptions and owns a registered
 * DG scalar field \c muAV.  The coefficient is cellwise constant:
 *
 * \f[
 *   \mu_{AV} = \chi(\rho)\,\mu_0\,\max_K(|U| + a)\,h/p,
 * \f]
 *
 * where \f$\chi\f$ is a modal smoothness activation based on the density field.
 *
 * Example \c constant/dgOptions syntax:
 *
 * \verbatim
 * artificialDiffusion
 * {
 *     enabled       true;
 *     type          NonSmooth;
 *     mu0           1.0;
 *     Skappa       -1.0;
 *     Kappa         0.25;
 *     sensorOffset  1;
 *     diffusiveCFL  0.25;
 *     report        false;
 * }
 * \endverbatim
 *
 * Dictionary entries:
 * - \c enabled: Switches the model on or off.  When false or when the
 *   \c artificialDiffusion block is absent, \c muAV is reset to zero.
 * - \c type: Artificial-diffusion model name.  The only supported value is
 *   \c NonSmooth.
 * - \c mu0: Multiplicative scaling for the artificial viscosity magnitude.
 *   Larger values add stronger artificial diffusion in troubled cells.
 * - \c Skappa: Base threshold for the base-10 logarithmic density sensor.
 *   The effective threshold is shifted by polynomial order as
 *   \f$S_\kappa - 4.25\log_{10}(p)\f$.
 * - \c Kappa: Half-width of the sine activation ramp.  The activation changes
 *   smoothly from 0 to 1 across
 *   \f$[S_\kappa-\kappa, S_\kappa+\kappa]\f$.
 * - \c sensorOffset: Number of highest modal shells included in the density
 *   smoothness sensor.  A value of 1 uses the top shell only.
 * - \c diffusiveCFL: CFL coefficient used to estimate the explicit diffusive
 *   time-step restriction from \c muAV.
 * - \c report: Prints per-cell sensor, activation, and \c muAV values for
 *   active cells when true.
 */
class dgArtificialDiffusionNonSmooth
:
    public IOdictionary
{
    // Private data

        /** \brief DG mesh supplying geometry, basis data, and p-order. */
        const dgGeomMesh& mesh_;

        /** \brief Density field used by the modal smoothness sensor. */
        const dgField<scalar>& rho_;

        /** \brief Velocity field used to estimate the local wave speed. */
        const dgField<vector>& U_;

        /** \brief Speed-of-sound field used to estimate the local wave speed. */
        const dgField<scalar>& a_;

        /** \brief Switch read from \c artificialDiffusion/enabled. */
        bool enabled_;

        /** \brief Artificial-diffusion model type; currently only \c NonSmooth. */
        word type_;

        /** \brief Registered DG field name for the artificial viscosity. */
        word muAVName_;

        /** \brief User scaling coefficient in the viscosity formula. */
        scalar mu0_;

        /** \brief Base logarithmic sensor threshold. */
        scalar Skappa_;

        /** \brief Half-width of the smooth sine activation interval. */
        scalar Kappa_;

        /** \brief Number of top modal shells included in the smoothness sensor. */
        label sensorOffset_;

        /** \brief CFL coefficient for artificial-diffusion time-step limits. */
        scalar diffusiveCFL_;

        /** \brief Emit per-cell sensor and viscosity diagnostics when active. */
        bool report_;

        /** \brief Owned artificial-viscosity DG scalar field. */
        dgField<scalar> muAV_;

        /** \brief Cached modal count for tetrahedral cells. */
        label tetNDof_;

        /** \brief Cached modal count for hexahedral cells. */
        label hexNDof_;

        /** \brief Cached modal count for prism cells. */
        label prismNDof_;

        /** \brief Cached modal count for pyramid cells. */
        label pyramidNDof_;

        /** \brief High-order modal indices for tetrahedral cells. */
        List<label> tetHighModeDof_;

        /** \brief High-order modal indices for hexahedral cells. */
        List<label> hexHighModeDof_;

        /** \brief High-order modal indices for prism cells. */
        List<label> prismHighModeDof_;

        /** \brief High-order modal indices for pyramid cells. */
        List<label> pyramidHighModeDof_;

        /**
         * \brief Build the IOobject used to read \c constant/dgOptions.
         *
         * The returned object uses \c MUST_READ_IF_MODIFIED when the file is
         * present, and \c NO_READ otherwise so the model can remain disabled
         * when no \c dgOptions file or \c artificialDiffusion block exists.
         *
         * \param mesh DG mesh whose fvMesh registry owns the dictionary.
         * \return IOobject configured for the dgOptions dictionary.
         */
        IOobject createIOobject(const dgGeomMesh& mesh) const;

        /**
         * \brief Return cached high-order modal indices for a cell topology.
         *
         * \param cellType Cell topology.
         * \return Cached list of high-order modal indices.
         */
        const List<label>& highModeDof(const dgCellType cellType) const;

        /**
         * \brief Return cached number of modal DoFs for a cell topology.
         *
         * \param cellType Cell topology.
         * \return Expected modal DoF count for the current p-order.
         */
        label expectedNDof(const dgCellType cellType) const;

        /**
         * \brief Compute the logarithmic density smoothness sensor for a cell.
         *
         * The sensor compares the L2 energy in the selected high-order modal
         * content against the L2 energy of the full density field.
         *
         * \param cellID Local cell index.
         * \return Base-10 logarithm of the normalized high-order content.
         */
        scalar densitySensor(const label cellID) const;

        /**
         * \brief Convert a logarithmic sensor value into an activation factor.
         *
         * Uses a sine ramp from 0 to 1 across
         * \f$[S_\kappa-\kappa, S_\kappa+\kappa]\f$.
         *
         * \param logSensor Base-10 logarithmic smoothness sensor.
         * \param pOrder Polynomial order used for p-dependent thresholding.
         * \return Activation coefficient in the range [0, 1].
         */
        scalar activation(const scalar logSensor, const label pOrder) const;

        /**
         * \brief Return the local length scale divided by polynomial order.
         *
         * \param cellID Local cell index.
         * \return \f$h/p\f$ scale used by the artificial viscosity formula.
         */
        scalar hOverP(const label cellID) const;

        /**
         * \brief Estimate maximum local acoustic wave speed in a cell.
         *
         * \param cellID Local cell index.
         * \return Maximum \f$|U| + a\f$ over cell Gauss points.
         */
        scalar maxWaveSpeed(const label cellID) const;

        /**
         * \brief Store a cellwise-constant artificial viscosity value.
         *
         * The cell mean mode is assigned \p mu and all higher modes are reset
         * to zero before refreshing the cell DoF representation.
         *
         * \param cellID Local cell index.
         * \param mu Artificial viscosity value for the cell.
         */
        void setCellMuAV(const label cellID, const scalar mu);

        /**
         * \brief Read \c artificialDiffusion coefficients from \c dgOptions.
         *
         * Recognized entries are \c enabled, \c type, \c mu0, \c Skappa,
         * \c Kappa, \c sensorOffset, \c diffusiveCFL, and \c report.
         */
        void readCoeffs();

public:

        TypeName("dgArtificialDiffusionNonSmooth");

        /**
         * \brief Construct the NonSmooth artificial-diffusion model.
         *
         * The constructor reads \c constant/dgOptions, creates the registered
         * \c muAV DG field, precomputes modal sensor index caches, and performs
         * an initial coefficient update.
         *
         * \param mesh DG mesh.
         * \param rho Density DG field.
         * \param U Velocity DG field.
         * \param a Speed-of-sound DG field.
         */
        dgArtificialDiffusionNonSmooth
        (
            const dgGeomMesh& mesh,
            const dgField<scalar>& rho,
            const dgField<vector>& U,
            const dgField<scalar>& a
        );

        /** \brief Copy construction is disabled. */
        dgArtificialDiffusionNonSmooth
        (
            const dgArtificialDiffusionNonSmooth&
        ) = delete;

        /** \brief Copy assignment is disabled. */
        dgArtificialDiffusionNonSmooth& operator=
        (
            const dgArtificialDiffusionNonSmooth&
        ) = delete;

        /**
         * \brief Return whether the model is enabled.
         *
         * \return True when \c artificialDiffusion/enabled is true and the
         *         selected type is supported.
         */
        inline bool enabled() const
        {
            return enabled_;
        }

        /**
         * \brief Recompute and return the artificial-viscosity field.
         *
         * When the model is disabled, the field is reset to zero and still
         * returned so callers can keep a uniform API.
         *
         * \return Mutable reference to the owned \c muAV field.
         */
        dgField<scalar>& updateMuAV();

        /**
         * \brief Compatibility wrapper for updateMuAV().
         */
        inline void update()
        {
            updateMuAV();
        }

        /**
         * \brief Return the explicit diffusive time-step limit for a cell.
         *
         * Uses the current \c muAV field and the configured \c diffusiveCFL.
         * Returns \c GREAT when the model is disabled or \c muAV is zero.
         *
         * \param cellID Local cell index.
         * \return Artificial-diffusion time-step restriction.
         */
        scalar diffusiveDeltaT(const label cellID) const;

        /**
         * \brief Return mutable access to the owned artificial-viscosity field.
         *
         * \return Mutable reference to \c muAV.
         */
        inline dgField<scalar>& muAV()
        {
            return muAV_;
        }

        /**
         * \brief Return read-only access to the owned artificial-viscosity field.
         *
         * \return Const reference to \c muAV.
         */
        inline const dgField<scalar>& muAV() const
        {
            return muAV_;
        }

        /**
         * \brief Write a short model summary to a stream.
         *
         * \param os Output stream.
         */
        void writeInfo(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
