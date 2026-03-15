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
    Foam::dgRK4TimeScheme

Description
    Classical explicit four-stage Runge-Kutta time integration for DG fields.

SourceFiles
    dgRK4TimeScheme.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgRK4TimeScheme_H
#define Foam_dgRK4TimeScheme_H

#include "dgTimeScheme.H"
#include "RK4FieldState.H"
#include "PtrList.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dgRK4TimeScheme Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgRK4TimeScheme
 * \brief Classical explicit RK4 time integration for DG fields.
 *
 * \details
 * The scheme owns one \c RK4FieldState per registered field and advances all
 * fields using the standard four-stage Runge-Kutta method
 *
 * \f[
 * \begin{aligned}
 * \mathbf{U}^{(1)} &= \mathbf{U}^{n} + \frac{1}{2}\Delta t\,\mathbf{M}^{-1}\mathbf{R}_1, \\
 * \mathbf{U}^{(2)} &= \mathbf{U}^{n} + \frac{1}{2}\Delta t\,\mathbf{M}^{-1}\mathbf{R}_2, \\
 * \mathbf{U}^{(3)} &= \mathbf{U}^{n} + \Delta t\,\mathbf{M}^{-1}\mathbf{R}_3, \\
 * \mathbf{U}^{n+1}
 * &= \mathbf{U}^{n}
 * + \frac{\Delta t}{6}
 * \left(
 *     \mathbf{M}^{-1}\mathbf{R}_1
 *   + 2\mathbf{M}^{-1}\mathbf{R}_2
 *   + 2\mathbf{M}^{-1}\mathbf{R}_3
 *   + \mathbf{M}^{-1}\mathbf{R}_4
 * \right).
 * \end{aligned}
 * \f]
 */
class dgRK4TimeScheme
:
    public dgTimeScheme
{
private:

        /** Names of registered scalar DG fields. */
        List<word> scalarNames_;

        /** Per-field RK4 storage for registered scalar fields. */
        PtrList<RK4FieldState<scalar>> scalarStates_;

        /** Names of registered vector DG fields. */
        List<word> vectorNames_;

        /** Per-field RK4 storage for registered vector fields. */
        PtrList<RK4FieldState<vector>> vectorStates_;

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

        /**
         * \brief Validate an explicit RK4 stage index.
         *
         * \param stageI Zero-based stage index
         */
        void validateStageIndex(const label stageI) const;

        /** \brief Validate the currently active RK4 stage index. */
        void validateStageIndex() const;

public:

    TypeName("RK4");

        /**
         * \brief Construct the RK4 time scheme.
         *
         * \param name Scheme instance name
         * \param dict Scheme dictionary
         * \param mesh DG mesh shared by all registered fields
         */
        dgRK4TimeScheme
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

        /** Disallow default construction. */
        dgRK4TimeScheme() = delete;

        /** Default destructor. */
        virtual ~dgRK4TimeScheme() = default;

        /** \brief Read RK4-specific settings from a dictionary. */
        virtual void read(const dictionary& dict) override;

        /** \brief Register a scalar DG field with the RK4 scheme. */
        virtual void registerField(dgField<scalar>& U) override;

        /** \brief Register a vector DG field with the RK4 scheme. */
        virtual void registerField(dgField<vector>& U) override;

        /** \brief Return the number of RK4 stages, namely 4. */
        virtual label nStages() const override;

        /** \brief Return a registered scalar solution field by name. */
        virtual dgField<scalar>& scalarField(const word& fieldName) override;

        /** \brief Return a registered scalar solution field by name. */
        virtual const dgField<scalar>& scalarField(const word& fieldName) const override;

        /** \brief Return a registered vector solution field by name. */
        virtual dgField<vector>& vectorField(const word& fieldName) override;

        /** \brief Return a registered vector solution field by name. */
        virtual const dgField<vector>& vectorField(const word& fieldName) const override;

        /** \brief Return the scalar residual buffer of one RK4 stage. */
        virtual List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) override;

        /** \brief Return the scalar residual buffer of one RK4 stage. */
        virtual const List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) const override;

        /** \brief Return the vector residual buffer of one RK4 stage. */
        virtual List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) override;

        /** \brief Return the vector residual buffer of one RK4 stage. */
        virtual const List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) const override;

        /** \brief Save one cell as the beginning-of-step RK4 state. */
        virtual void saveOldFields(const label cellID) override;

        /** \brief Restore one cell from its beginning-of-step RK4 state. */
        virtual void restoreOldFields(const label cellID) override;

        /** \brief RK4 stores its stage residuals only, so this is a no-op. */
        virtual void saveStageFields(const label cellID) override;

        /** \brief RK4 stores its stage residuals only, so this is a no-op. */
        virtual void restoreStageFields(const label cellID) override;

        /** \brief Clear the residuals of one cell for all registered fields. */
        virtual void clearResiduals(const label cellID) override;

        /** \brief Clear the residuals of all cells for all registered fields. */
        virtual void clearResiduals() override;

        /**
         * \brief Apply the RK4 stage update to one cell.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        virtual void updateStage
        (
            const label cellID,
            const scalar dt
        ) override;

        /** \brief RK4 completes its weighted combination in the last stage, so no extra work is needed. */
        virtual void finalizeTimeStep
        (
            const label cellID,
            const scalar dt
        ) override;

        /** \brief Reset all RK4 field states. */
        virtual void reset() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
