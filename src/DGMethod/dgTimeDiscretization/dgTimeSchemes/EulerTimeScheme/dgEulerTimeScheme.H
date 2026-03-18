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
    Foam::dgEulerTimeScheme

Description
    Explicit one-stage Euler time integration for DG fields.

SourceFiles
    dgEulerTimeScheme.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgEulerTimeScheme_H
#define Foam_dgEulerTimeScheme_H

#include "dgTimeScheme.H"
#include "EulerFieldState.H"
#include "PtrList.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dgEulerTimeScheme Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgEulerTimeScheme
 * \brief Explicit one-stage Euler time integration for DG fields.
 *
 * \details
 * The scheme owns one \c EulerFieldState per registered field and advances all
 * fields using
 *
 * \f[
 *   \mathbf{U}^{n+1}
 *   =
 *   \mathbf{U}^{n}
 *   +
 *   \Delta t\,\mathbf{M}^{-1}\mathbf{R}^{n}.
 * \f]
 *
 * Since Euler has only one stage, the active stage index must always be zero
 * and no intermediate-stage field storage is required.
 */
class dgEulerTimeScheme
:
    public dgTimeScheme
{
private:

        /** Names of registered scalar DG fields. */
        List<word> scalarNames_;

        /** Per-field Euler storage for registered scalar fields. */
        PtrList<EulerFieldState<scalar>> scalarStates_;

        /** Names of registered vector DG fields. */
        List<word> vectorNames_;

        /** Per-field Euler storage for registered vector fields. */
        PtrList<EulerFieldState<vector>> vectorStates_;

        /**
         * \brief Return the index of a registered scalar field.
         *
         * \param fieldName Registered field name
         * \return State index, or \c -1 if not found
         */
        label findScalarIndex(const word& fieldName) const;

        /**
         * \brief Return the index of a registered vector field.
         *
         * \param fieldName Registered field name
         * \return State index, or \c -1 if not found
         */
        label findVectorIndex(const word& fieldName) const;

        /** \brief Check that Euler is used only with stage index zero. */
        void validateStageIndex() const;

public:

    TypeName("Euler");

        /**
         * \brief Construct the Euler time scheme.
         *
         * \param name Scheme instance name
         * \param dict Scheme dictionary
         * \param mesh DG mesh shared by all registered fields
         */
        dgEulerTimeScheme
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

        /** Disallow default construction. */
        dgEulerTimeScheme() = delete;

        /** Default destructor. */
        virtual ~dgEulerTimeScheme() = default;

        /** \brief Read Euler-specific settings from a dictionary. */
        virtual void read(const dictionary& dict) override;

        /** \brief Register a scalar DG field with the Euler scheme. */
        virtual void registerField(dgField<scalar>& U) override;

        /** \brief Register a vector DG field with the Euler scheme. */
        virtual void registerField(dgField<vector>& U) override;

        /** \brief Return the number of stages of explicit Euler, namely 1. */
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

        /** \brief Return the scalar residual buffer of the current Euler stage. */
        virtual List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) override;

        /** \brief Return the scalar residual buffer of the current Euler stage. */
        virtual const List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) const override;

        /** \brief Return the vector residual buffer of the current Euler stage. */
        virtual List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) override;

        /** \brief Return the vector residual buffer of the current Euler stage. */
        virtual const List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) const override;

        /** \brief Save one cell as the beginning-of-step Euler state. */
        virtual void saveOldFields(const label cellID) override;

        /** \brief Restore one cell from its beginning-of-step Euler state. */
        virtual void restoreOldFields(const label cellID) override;

        /** \brief Euler has no persistent intermediate-stage fields, so this is a no-op. */
        virtual void saveStageFields(const label cellID) override;

        /** \brief Euler has no persistent intermediate-stage fields, so this is a no-op. */
        virtual void restoreStageFields(const label cellID) override;

        /** \brief Clear the residuals of one cell for all registered fields. */
        virtual void clearResiduals(const label cellID) override;

        /** \brief Clear the residuals of all cells for all registered fields. */
        virtual void clearResiduals() override;

        /**
         * \brief Apply the explicit Euler update to one cell.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        virtual void updateStage
        (
            const label cellID,
            const scalar dt
        ) override;

        /** \brief Euler completes its update in the only stage, so no extra work is needed. */
        virtual void finalizeTimeStep
        (
            const label cellID,
            const scalar dt
        ) override;

        /** \brief Reset all Euler field states. */
        virtual void reset() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
