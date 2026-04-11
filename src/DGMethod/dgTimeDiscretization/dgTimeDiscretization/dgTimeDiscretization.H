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
    Foam::dgTimeDiscretization

Description
    Manager for DG time discretization selected from the dgSchemes dictionary.

SourceFiles
    dgTimeDiscretization.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgTimeDiscretization_H
#define Foam_dgTimeDiscretization_H

#include "Time.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "dgGeomMesh.H"
#include "dgGeneralPDETerm.H"
#include "dgTimeScheme.H"
#include "error.H"

#include <type_traits>

namespace Foam
{

class Ostream;

/*---------------------------------------------------------------------------*\
                    Class dgTimeDiscretization Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgTimeDiscretization
 * \brief High-level driver that manages the explicit stage loop for DG fields.
 *
 * \details
 * This class owns the runtime-selected concrete \c dgTimeScheme and exposes a
 * uniform driver interface to the solver. The stage loop is advanced through
 * \c loop(), residuals are accumulated through \c addResidual(), and the
 * per-cell explicit update is applied through \c updateStage().
 *
 * At the semi-discrete level, the solver advances
 *
 * \f[
 *     \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{R}(\mathbf{U}),
 * \f]
 *
 * where the concrete scheme determines how the stage states
 * \f$\mathbf{U}^{(s)}\f$ are combined. The present driver is agnostic to the
 * actual stage formula and delegates that work to the selected scheme.
 *
 * Operationally, one physical time step follows
 *
 * \f[
 *     \mathbf{U}^n
 *     \xrightarrow{\text{save old state}}
 *     \mathbf{U}^{(0)}, \mathbf{U}^{(1)}, \dots
 *     \xrightarrow{\text{final stage}}
 *     \mathbf{U}^{n+1}.
 * \f]
 */
class dgTimeDiscretization
{
private:

        /** Helper used to trigger a clean compile-time error for unsupported residual types. */
        template<class Type>
        struct alwaysFalse : std::false_type
        {};

        /**
         * \brief Return the residual buffer of the current stage for a given field type.
         *
         * \param fieldName Registered DG field name
         * \return Residual storage indexed by [cell][dof]
         */
        template<class Type>
        inline List<List<Type>>& residualBuffer(const word& fieldName)
        {
            if constexpr (std::is_same_v<Type, scalar>)
            {
                return scalarResidual(fieldName, stage());
            }
            else if constexpr (std::is_same_v<Type, vector>)
            {
                return vectorResidual(fieldName, stage());
            }
            else
            {
                static_assert(alwaysFalse<Type>::value, "Unsupported residual type");
            }
        }

        /** OpenFOAM time database providing \f$\Delta t\f$ and write control. */
        Time& runTime_;

        /** Shared DG mesh used when iterating over all cells. */
        const dgGeomMesh& mesh_;

        /** Cached \c ddtSchemes sub-dictionary extracted from \c dgSchemes. */
        dictionary ddtSchemeDict_;

        /** Name of the selected explicit time scheme, for example \c Euler or \c SSPRK3. */
        word schemeType_;

        /** Runtime-selected concrete time-scheme object. */
        autoPtr<dgTimeScheme> timeScheme_;

        /** Total number of stages of the selected scheme. */
        label nStage_;

        /** Current zero-based stage index seen by the driver. */
        label stageI_;

        /** True once the current physical time step has entered its first stage. */
        bool timeStepStarted_;

        /** Use per-cell pseudo-time increments instead of the global time-step. */
        Switch localTimeStepping_;

        /** Maximum allowed ratio between local pseudo-time and global minimum time-step. */
        scalar localTimeStepMaxRatio_;

        /** Per-cell pseudo-time increments used when local time stepping is enabled. */
        List<scalar> localDeltaT_;

        /**
         * Reference scalar residual summaries built from the warmup time steps.
         *
         * Indexed as [field][stage][dof]. Empty stage entries mean that the
         * corresponding reference has not been stored yet.
         */
        mutable List<List<List<scalar>>> scalarResidualReference_;

        /**
         * Scalar reference scales for vector residual summaries.
         *
         * Indexed as [field][stage][dof]. Each DoF stores the maximum over the
         * three vector components and the warmup time steps, so later vector
         * residuals are scaled by one shared scalar reference.
         */
        mutable List<List<List<scalar>>> vectorResidualReference_;

        /** Number of physical time steps used to build the residual reference. */
        label residualReferenceWarmupSteps_;

        /** Number of physical time steps completed since startup. */
        label completedTimeSteps_;

public:

    // Constructors

        /**
         * \brief Construct the DG time-discretization driver.
         *
         * The manager expects a dictionary block of the form
         * \verbatim
         * ddtSchemes
         * {
         *     scheme SSPRK3;
         * }
         * \endverbatim
         *
         * \param dgSchemesDict Top-level DG discretization dictionary
         * \param runTime OpenFOAM time database
         * \param mesh DG mesh used by the solver
         */
        dgTimeDiscretization
        (
            const dictionary& dgSchemesDict,
            Time& runTime,
            const dgGeomMesh& mesh
        );

        /** Disallow default construction. */
        dgTimeDiscretization() = delete;

        /** Default destructor. */
        ~dgTimeDiscretization() = default;

    // Access

        /** \brief Return the selected time-scheme type name. */
        inline const word& schemeType() const
        {
            return schemeType_;
        }

        /** \brief Return read-only access to the OpenFOAM run-time object. */
        inline const Time& runTime() const
        {
            return runTime_;
        }

        /** \brief Return writable access to the OpenFOAM run-time object. */
        inline Time& runTime()
        {
            return runTime_;
        }

        /** \brief Return the current physical time-step as a dimensioned scalar. */
        inline dimensionedScalar deltaT() const
        {
            return runTime_.deltaT();
        }

        /** \brief Return the current physical time-step value. */
        inline scalar deltaTValue() const
        {
            return runTime_.deltaTValue();
        }

        /** \brief Return true when pseudo-local time stepping is enabled. */
        inline bool localTimeStepping() const
        {
            return localTimeStepping_;
        }

        /** \brief Return the maximum local/global pseudo-time-step ratio. */
        inline scalar localTimeStepMaxRatio() const
        {
            return localTimeStepMaxRatio_;
        }

        /**
         * \brief Return the time increment used to update one cell.
         *
         * With local pseudo-time stepping disabled, this is the global
         * \c runTime.deltaTValue(). With LTS enabled, this is the capped
         * per-cell pseudo-time increment.
         */
        inline scalar cellDeltaT(const label cellID) const
        {
            if (!localTimeStepping_)
            {
                return deltaTValue();
            }

            if (cellID < 0 || cellID >= localDeltaT_.size())
            {
                FatalErrorInFunction
                    << "Cell index " << cellID
                    << " is out of range for localDeltaT of size "
                    << localDeltaT_.size() << '.'
                    << abort(FatalError);
            }

            return localDeltaT_[cellID];
        }

        /** \brief Return the cached \c ddtSchemes dictionary block. */
        inline const dictionary& dict() const
        {
            return ddtSchemeDict_;
        }

        /** \brief Return the number of explicit stages of the active time scheme. */
        inline label nStages() const
        {
            return nStage_;
        }

        /** \brief Alias of \c nStages(). */
        inline label nStage() const
        {
            return nStage_;
        }

        /** \brief Return read-only access to the active time-scheme object. */
        inline const dgTimeScheme& scheme() const
        {
            return timeScheme_();
        }

        /** \brief Return writable access to the active time-scheme object. */
        inline dgTimeScheme& scheme()
        {
            return timeScheme_();
        }

    // Registration

        /**
         * \brief Register a scalar DG field with the active time scheme.
         *
         * \param U Scalar DG solution field
         */
        inline void registerField(dgField<scalar>& U)
        {
            timeScheme_().registerField(U);
        }

        /**
         * \brief Register a vector DG field with the active time scheme.
         *
         * \param U Vector DG solution field
         */
        inline void registerField(dgField<vector>& U)
        {
            timeScheme_().registerField(U);
        }

    // Driver API

        /**
         * \brief Set a single global explicit time step.
         *
         * This updates the OpenFOAM runtime clock step. If local pseudo-time
         * stepping is enabled, all local increments are also reset to this
         * global value.
         */
        void setDeltaT(const scalar dt);

        /**
         * \brief Set candidate cell-wise explicit time steps.
         *
         * The OpenFOAM runtime clock is set to the global minimum of
         * \p cellDeltaT. If local pseudo-time stepping is enabled, the
         * selected time scheme uses the capped per-cell values while the
         * runtime clock still advances by the global minimum.
         */
        void setDeltaT(const List<scalar>& cellDeltaT);

        /**
         * \brief Set the current stage index on both the driver and the active scheme.
         *
         * \param stageI Zero-based stage index
         */
        inline void setStage(const label stageI)
        {
            stageI_ = stageI;
            timeScheme_().setStage(stageI);
        }

        /** \brief Return the current stage index managed by the driver. */
        inline label stage() const
        {
            return stageI_;
        }

        /**
         * \brief Advance the stage-loop state machine.
         *
         * The first call of a physical time step initializes stage 0 and saves
         * the old solution fields. Later calls advance to the next stage until
         * all stages are completed.
         *
         * \return \c true while there is a stage to execute, otherwise \c false
         */
        inline bool loop()
        {
            if (!timeStepStarted_)
            {
                setStage(0);
                saveOldFields();

                timeStepStarted_ = true;
                return true;
            }

            if (stageI_ + 1 >= nStage_)
            {
                return false;
            }

            setStage(stageI_ + 1);
            return true;
        }

        /** \brief Save the old solution of all cells at the start of a physical time step. */
        inline void saveOldFields()
        {
            for (label cellID = 0; cellID < mesh_.nCells(); ++cellID)
            {
                timeScheme_().saveOldFields(cellID);
            }
        }

        /**
         * \brief Restore the previously saved solution of one cell.
         *
         * \param cellID Local DG cell index
         */
        inline void restoreOldFields(const label cellID)
        {
            timeScheme_().restoreOldFields(cellID);
        }

        /**
         * \brief Save the current intermediate stage state of one cell.
         *
         * \param cellID Local DG cell index
         */
        inline void saveStageFields(const label cellID)
        {
            timeScheme_().saveStageFields(cellID);
        }

        /**
         * \brief Restore a previously saved intermediate stage state of one cell.
         *
         * \param cellID Local DG cell index
         */
        inline void restoreStageFields(const label cellID)
        {
            timeScheme_().restoreStageFields(cellID);
        }

        /**
         * \brief Clear residual storage of one cell in the active scheme.
         *
         * \param cellID Local DG cell index
         */
        inline void clearResiduals(const label cellID)
        {
            timeScheme_().clearResiduals(cellID);
        }

        /** \brief Clear all residual storage of the active scheme. */
        inline void clearResiduals()
        {
            timeScheme_().clearResiduals();
        }

        /**
         * \brief Apply one stage update to a single cell with an explicit time-step.
         *
         * On the last stage, this also finalizes the physical time step of the
         * given cell.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        inline void updateStage(const label cellID, const scalar dt)
        {
            timeScheme_().updateStage(cellID, dt);

            if (stageI_ + 1 == nStage_)
            {
                timeScheme_().finalizeTimeStep(cellID, dt);
            }
        }

        /**
         * \brief Apply one stage update to a single cell using the active time increment.
         *
         * \param cellID Local DG cell index
         */
        inline void updateStage(const label cellID)
        {
            updateStage(cellID, cellDeltaT(cellID));
        }

        /**
         * \brief Finalize the physical time step of one cell with an explicit \p dt.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        inline void finalizeTimeStep(const label cellID, const scalar dt)
        {
            timeScheme_().finalizeTimeStep(cellID, dt);
        }

        /**
         * \brief Finalize the physical time step of one cell using \c runTime.deltaTValue().
         *
         * \param cellID Local DG cell index
         */
        inline void finalizeTimeStep(const label cellID)
        {
            timeScheme_().finalizeTimeStep(cellID, cellDeltaT(cellID));
        }

        /**
         * \brief Reset the driver after a physical time step has completed.
         *
         * This clears accumulated residuals and marks the next \c loop() call as
         * the beginning of a new physical time step.
         */
        inline void reset()
        {
            setStage(0);
            clearResiduals();
            timeStepStarted_ = false;
            ++completedTimeSteps_;
        }

    // Field/residual access

        /**
         * \brief Return a registered scalar DG field by name.
         *
         * \param fieldName Registered field name
         * \return Writable scalar field reference
         */
        inline dgField<scalar>& scalarField(const word& fieldName)
        {
            return timeScheme_().scalarField(fieldName);
        }

        /**
         * \brief Return a registered scalar DG field by name.
         *
         * \param fieldName Registered field name
         * \return Read-only scalar field reference
         */
        inline const dgField<scalar>& scalarField(const word& fieldName) const
        {
            return timeScheme_().scalarField(fieldName);
        }

        /**
         * \brief Return a registered vector DG field by name.
         *
         * \param fieldName Registered field name
         * \return Writable vector field reference
         */
        inline dgField<vector>& vectorField(const word& fieldName)
        {
            return timeScheme_().vectorField(fieldName);
        }

        /**
         * \brief Return a registered vector DG field by name.
         *
         * \param fieldName Registered field name
         * \return Read-only vector field reference
         */
        inline const dgField<vector>& vectorField(const word& fieldName) const
        {
            return timeScheme_().vectorField(fieldName);
        }

        /**
         * \brief Return the scalar residual buffer for a given field and stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         * \return Writable scalar residual storage
         */
        inline List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        )
        {
            return timeScheme_().scalarResidual(fieldName, stageI);
        }

        /**
         * \brief Return the scalar residual buffer for a given field and stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         * \return Read-only scalar residual storage
         */
        inline const List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) const
        {
            return timeScheme_().scalarResidual(fieldName, stageI);
        }

        /**
         * \brief Return the vector residual buffer for a given field and stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         * \return Writable vector residual storage
         */
        inline List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        )
        {
            return timeScheme_().vectorResidual(fieldName, stageI);
        }

        /**
         * \brief Return the vector residual buffer for a given field and stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         * \return Read-only vector residual storage
         */
        inline const List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) const
        {
            return timeScheme_().vectorResidual(fieldName, stageI);
        }

        /**
         * \brief Accumulate one residual contribution into the current-stage buffer of a cell.
         *
         * \param fieldName Registered field name
         * \param cellID Local DG cell index
         * \param residual Residual vector to accumulate
         */
        template<class Type>
        inline void addResidual
        (
            const word& fieldName,
            const label cellID,
            const List<Type>& residual
        )
        {
            List<Type>& cellResidual = residualBuffer<Type>(fieldName)[cellID];

            forAll(cellResidual, dofI)
            {
                cellResidual[dofI] += residual[dofI];
            }
        }

        /**
         * \brief Write a short summary of the active time-discretization settings.
         *
         * \param os Output stream
         */
        void writeInfo(Ostream& os) const;

        /**
         * \brief Write optional time-discretization parameters.
         *
         * Parameters that are inactive or not meaningful are intentionally not
         * printed. For example, local pseudo-time-stepping options are printed
         * only when local time stepping is enabled.
         *
         * \param os Output stream
         */
        void writeParameters(Ostream& os) const;

        /**
         * \brief Write the normalized maximum absolute residual of every DOF.
         *
         * The first \c residualReferenceWarmupSteps_ physical time steps are
         * used to build a per-field, per-stage, per-DoF maximum reference.
         * Later outputs are divided by that warmup reference.
         *
         * \param os Output stream
         */
        void writeResiduals(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
