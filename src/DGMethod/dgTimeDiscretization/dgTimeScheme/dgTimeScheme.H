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
    Foam::dgTimeScheme

Description
    Base wrapper for DG time-integration schemes.

    dgTimeScheme provides the runtime-selectable interface shared by
    explicit time-integration schemes such as Euler, RK4, SSPRK2, SSPRK3,
    and SSPRK4. Derived classes own scheme-specific field-state types and
    implement the actual stage-update formulas.

SourceFiles
    dgTimeScheme.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgTimeScheme_H
#define Foam_dgTimeScheme_H

#include "word.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "scalar.H"
#include "vector.H"
#include "List.H"
#include "dgGeomMesh.H"
#include "dgField.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class dgTimeScheme Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgTimeScheme
 * \brief Abstract interface for explicit DG time-integration schemes.
 *
 * \details
 * A concrete time scheme advances the semi-discrete DG system
 *
 * \f[
 *     \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{R}(\mathbf{U}),
 * \f]
 *
 * by applying a scheme-specific explicit stage update to each registered DG
 * field. In the current solver design, the application assembles stage
 * residuals per cell and the concrete scheme maps them into the next stage
 * through formulas of the generic form
 *
 * \f[
 *     \mathbf{U}^{(s+1)} = \mathcal{F}_s
 *     \left(
 *         \mathbf{U}^n,\,
 *         \mathbf{U}^{(0)},\dots,\mathbf{U}^{(s)},\,
 *         \Delta t,\,
 *         \mathbf{M}^{-1}\mathbf{R}^{(s)}
 *     \right).
 * \f]
 *
 * The base class itself does not store residuals or intermediate stage
 * values. Those responsibilities are delegated to the derived classes through
 * their field-state containers.
 */
class dgTimeScheme
{
protected:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Name used to register this scheme instance in the runtime-selection framework. */
        word name_;

        /** Local dictionary containing the coefficients or options of the concrete scheme. */
        dictionary dict_;

        /** Shared DG mesh used to access cell-local geometric data and mass-matrix diagonals. */
        const dgGeomMesh& mesh_;

        /** Zero-based index of the currently active explicit stage. */
        label stageI_;

public:

    //====================================================================//
    //                     Runtime type information
    //====================================================================//

    TypeName("dgTimeScheme");

    declareRunTimeSelectionTable
    (
        autoPtr,
        dgTimeScheme,
        dictionary,
        (const word& name, const dictionary& dict, const dgGeomMesh& mesh),
        (name, dict, mesh)
    );


    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct a time-scheme wrapper.
         *
         * \param name Scheme instance name
         * \param dict Dictionary containing scheme settings
         * \param mesh DG mesh shared by all registered fields
         */
        dgTimeScheme
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

        /** Disallow default construction. */
        dgTimeScheme() = delete;

        /** Virtual destructor. */
        virtual ~dgTimeScheme() = default;


    //====================================================================//
    //                              Factory
    //====================================================================//

        /**
         * \brief Select and construct a concrete time scheme.
         *
         * \param name Scheme instance name
         * \param timeSchemeType Selected time-scheme type name
         * \param dict Scheme dictionary
         * \param mesh DG geometric mesh
         * \return Runtime-selected concrete scheme
         */
        static autoPtr<dgTimeScheme> New
        (
            const word& name,
            const word& timeSchemeType,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );


    //====================================================================//
    //                              Access
    //====================================================================//

        /** \brief Return the runtime-selected name of this scheme instance. */
        inline const word& name() const
        {
            return name_;
        }

        /** \brief Return the local coefficient dictionary of this scheme. */
        inline const dictionary& dict() const
        {
            return dict_;
        }

        /** \brief Return the DG mesh associated with this scheme. */
        inline const dgGeomMesh& mesh() const
        {
            return mesh_;
        }

        /**
         * \brief Set the currently active explicit stage index.
         *
         * \param stageI Zero-based stage index
         */
        inline void setStage(const label stageI)
        {
            stageI_ = stageI;
        }

        /** \brief Return the zero-based index of the current explicit stage. */
        inline label stage() const
        {
            return stageI_;
        }


    //====================================================================//
    //                              Setup
    //====================================================================//

        /**
         * \brief Read or reread scheme settings from a dictionary.
         *
         * \param dict Scheme-specific dictionary block
         */
        virtual void read(const dictionary& dict) = 0;

        /**
         * \brief Register a scalar DG solution field with the scheme.
         *
         * \param U Scalar DG field that will be advanced in time
         */
        virtual void registerField(dgField<scalar>& U) = 0;

        /**
         * \brief Register a vector DG solution field with the scheme.
         *
         * \param U Vector DG field that will be advanced in time
         */
        virtual void registerField(dgField<vector>& U) = 0;


    //====================================================================//
    //                              Accessors
    //====================================================================//

        /** \brief Return the total number of explicit stages of the scheme. */
        virtual label nStages() const = 0;

        /**
         * \brief Return a registered scalar solution field by name.
         *
         * \param fieldName Registered DG field name
         * \return Writable scalar field reference
         */
        virtual dgField<scalar>& scalarField(const word& fieldName) = 0;

        /**
         * \brief Return a registered scalar solution field by name.
         *
         * \param fieldName Registered DG field name
         * \return Read-only scalar field reference
         */
        virtual const dgField<scalar>& scalarField(const word& fieldName) const = 0;

        /**
         * \brief Return a registered vector solution field by name.
         *
         * \param fieldName Registered DG field name
         * \return Writable vector field reference
         */
        virtual dgField<vector>& vectorField(const word& fieldName) = 0;

        /**
         * \brief Return a registered vector solution field by name.
         *
         * \param fieldName Registered DG field name
         * \return Read-only vector field reference
         */
        virtual const dgField<vector>& vectorField(const word& fieldName) const = 0;

        /**
         * \brief Return one scalar residual buffer for a given stage.
         *
         * \param fieldName Registered DG field name
         * \param stageI Zero-based stage index
         * \return Writable residual buffer indexed by [cell][dof]
         */
        virtual List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) = 0;

        /**
         * \brief Return one scalar residual buffer for a given stage.
         *
         * \param fieldName Registered DG field name
         * \param stageI Zero-based stage index
         * \return Read-only residual buffer indexed by [cell][dof]
         */
        virtual const List<List<scalar>>& scalarResidual
        (
            const word& fieldName,
            const label stageI
        ) const = 0;

        /**
         * \brief Return one vector residual buffer for a given stage.
         *
         * \param fieldName Registered DG field name
         * \param stageI Zero-based stage index
         * \return Writable residual buffer indexed by [cell][dof]
         */
        virtual List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) = 0;

        /**
         * \brief Return one vector residual buffer for a given stage.
         *
         * \param fieldName Registered DG field name
         * \param stageI Zero-based stage index
         * \return Read-only residual buffer indexed by [cell][dof]
         */
        virtual const List<List<vector>>& vectorResidual
        (
            const word& fieldName,
            const label stageI
        ) const = 0;


    //====================================================================//
    //                              Core API
    //====================================================================//

        /**
         * \brief Save the current solution of one cell as the beginning-of-step state.
         *
         * \param cellID Local DG cell index
         */
        virtual void saveOldFields(const label cellID) = 0;

        /**
         * \brief Restore one cell from its beginning-of-step state.
         *
         * \param cellID Local DG cell index
         */
        virtual void restoreOldFields(const label cellID) = 0;

        /**
         * \brief Save the current intermediate stage state of one cell, if required.
         *
         * Schemes without persistent intermediate-stage storage may implement
         * this as a no-op.
         *
         * \param cellID Local DG cell index
         */
        virtual void saveStageFields(const label cellID) = 0;

        /**
         * \brief Restore one previously saved intermediate stage state, if available.
         *
         * Schemes without persistent intermediate-stage storage may implement
         * this as a no-op.
         *
         * \param cellID Local DG cell index
         */
        virtual void restoreStageFields(const label cellID) = 0;

        /**
         * \brief Clear all stage residuals of one cell for all registered fields.
         *
         * \param cellID Local DG cell index
         */
        virtual void clearResiduals(const label cellID) = 0;

        /** \brief Clear all residual buffers over all cells, stages, and registered fields. */
        virtual void clearResiduals() = 0;

        /**
         * \brief Apply the scheme-specific update for one cell at the current stage.
         *
         * The residual buffer of the current stage is assumed to already contain
         * the assembled semi-discrete DG operator for \p cellID.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        virtual void updateStage
        (
            const label cellID,
            const scalar dt
        ) = 0;

        /**
         * \brief Finalize one cell after the last explicit stage.
         *
         * \param cellID Local DG cell index
         * \param dt Physical time-step size
         */
        virtual void finalizeTimeStep
        (
            const label cellID,
            const scalar dt
        ) = 0;

        /** \brief Reset the internal field-state storage of the scheme. */
        virtual void reset() = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
