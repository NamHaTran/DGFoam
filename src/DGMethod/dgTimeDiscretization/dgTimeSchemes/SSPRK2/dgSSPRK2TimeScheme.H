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
    Foam::dgSSPRK2TimeScheme

Description
    Two-stage second-order strong-stability-preserving Runge-Kutta time
    integration for DG fields.

SourceFiles
    dgSSPRK2TimeScheme.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgSSPRK2TimeScheme_H
#define Foam_dgSSPRK2TimeScheme_H

#include "dgTimeScheme.H"
#include "SSPRK2FieldState.H"
#include "PtrList.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgSSPRK2TimeScheme Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgSSPRK2TimeScheme
 * \brief Explicit SSPRK2 time scheme for DG fields.
 *
 * \details
 * The scheme advances the semi-discrete DG system
 *
 * \f[
 *     \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{R}(\mathbf{U})
 * \f]
 *
 * with the Shu-Osher SSPRK2 update
 *
 * \f[
 * \begin{aligned}
 * \mathbf{U}^{(1)} &=
 *     \mathbf{U}^{n} + \Delta t\,\mathbf{M}^{-1}\mathbf{R}_1, \\
 * \mathbf{U}^{n+1} &=
 *     \frac{1}{2}\mathbf{U}^{n}
 *   + \frac{1}{2}
 *     \left(
 *         \mathbf{U}^{(1)}
 *       + \Delta t\,\mathbf{M}^{-1}\mathbf{R}_2
 *     \right).
 * \end{aligned}
 * \f]
 *
 * One \c SSPRK2FieldState is stored for each registered scalar and vector
 * field, together with two stage residual buffers. The driver computes the
 * current residual, and this class converts it into the next stage state
 * using the diagonal DG mass matrix on each cell.
 */
class dgSSPRK2TimeScheme
:
    public dgTimeScheme
{
private:

        /** Names of registered scalar DG fields. */
        List<word> scalarNames_;

        /** Per-field SSPRK2 storage for scalar DG fields. */
        PtrList<SSPRK2FieldState<scalar>> scalarStates_;

        /** Names of registered vector DG fields. */
        List<word> vectorNames_;

        /** Per-field SSPRK2 storage for vector DG fields. */
        PtrList<SSPRK2FieldState<vector>> vectorStates_;

        /**
         * \brief Return the storage index of a registered scalar field.
         *
         * \param fieldName Registered field name
         * \return State index, or \c -1 if the field is not registered
         */
        label findScalarIndex(const word& fieldName) const;

        /**
         * \brief Return the storage index of a registered vector field.
         *
         * \param fieldName Registered field name
         * \return State index, or \c -1 if the field is not registered
         */
        label findVectorIndex(const word& fieldName) const;

        /**
         * \brief Validate an explicit SSPRK2 stage index.
         *
         * \param stageI Zero-based stage index, expected in \c [0,1]
         */
        void validateStageIndex(const label stageI) const;

        /** \brief Validate the currently active stage index stored in the base class. */
        void validateStageIndex() const;

public:

    TypeName("SSPRK2");

        /**
         * \brief Construct the SSPRK2 time scheme.
         *
         * \param name Scheme instance name
         * \param dict Scheme dictionary
         * \param mesh DG mesh shared by all registered fields
         */
        dgSSPRK2TimeScheme
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

        /** Disallow default construction. */
        dgSSPRK2TimeScheme() = delete;

        /** Default destructor. */
        virtual ~dgSSPRK2TimeScheme() = default;

        /** \brief Read SSPRK2-specific settings from a dictionary. */
        virtual void read(const dictionary& dict) override;

        /** \brief Register a scalar DG field with the time scheme. */
        virtual void registerField(dgField<scalar>& U) override;

        /** \brief Register a vector DG field with the time scheme. */
        virtual void registerField(dgField<vector>& U) override;

        /** \brief Return the number of SSPRK2 stages, namely 2. */
        virtual label nStages() const override;

        /** \brief Return a registered scalar solution field by name. */
        virtual dgField<scalar>& scalarField(const word& fieldName) override;

        /** \brief Return a registered scalar solution field by name. */
        virtual const dgField<scalar>& scalarField(const word& fieldName) const override;

        /** \brief Return a registered vector solution field by name. */
        virtual dgField<vector>& vectorField(const word& fieldName) override;

        /** \brief Return a registered vector solution field by name. */
        virtual const dgField<vector>& vectorField(const word& fieldName) const override;

        /** \brief Return the names of all registered scalar DG fields. */
        virtual const List<word>& scalarFieldNames() const override
        {
            return scalarNames_;
        }

        /** \brief Return the names of all registered vector DG fields. */
        virtual const List<word>& vectorFieldNames() const override
        {
            return vectorNames_;
        }

        /**
         * \brief Return the scalar residual buffer of one SSPRK2 stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         */
        virtual List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) override;

        /**
         * \brief Return the scalar residual buffer of one SSPRK2 stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         */
        virtual const List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) const override;

        /**
         * \brief Return the vector residual buffer of one SSPRK2 stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         */
        virtual List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) override;

        /**
         * \brief Return the vector residual buffer of one SSPRK2 stage.
         *
         * \param fieldName Registered field name
         * \param stageI Zero-based stage index
         */
        virtual const List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) const override;

        /**
         * \brief Save the beginning-of-step state for one cell.
         *
         * \param cellID Local DG cell index
         */
        virtual void saveOldFields(const label cellID) override;

        /**
         * \brief Restore the beginning-of-step state for one cell.
         *
         * \param cellID Local DG cell index
         */
        virtual void restoreOldFields(const label cellID) override;

        /**
         * \brief Save the current stage state for one cell.
         *
         * In SSPRK2 this stores \f$\mathbf{U}^{(1)}\f$ after the first stage.
         *
         * \param cellID Local DG cell index
         */
        virtual void saveStageFields(const label cellID) override;

        /**
         * \brief Restore the previously saved stage state for one cell.
         *
         * \param cellID Local DG cell index
         */
        virtual void restoreStageFields(const label cellID) override;

        /**
         * \brief Clear all stage residuals for one cell.
         *
         * \param cellID Local DG cell index
         */
        virtual void clearResiduals(const label cellID) override;

        /** \brief Clear all stage residuals on all cells for all registered fields. */
        virtual void clearResiduals() override;

        /**
         * \brief Apply the current SSPRK2 stage update to one cell.
         *
         * The method expects that the residual buffer for the active stage has
         * already been assembled by the solver, and then forms the next state
         * using \f$\mathbf{M}^{-1}\mathbf{R}\f$ with the cell-local diagonal
         * mass matrix.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        virtual void updateStage
        (
            const label cellID,
            const scalar dt
        ) override;

        /**
         * \brief Finalize the physical time step for one cell.
         *
         * SSPRK2 completes the final solution update in stage 1, so this is a
         * no-op hook required by the common interface.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        virtual void finalizeTimeStep
        (
            const label cellID,
            const scalar dt
        ) override;

        /** \brief Reset all stored field snapshots and residual buffers. */
        virtual void reset() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
