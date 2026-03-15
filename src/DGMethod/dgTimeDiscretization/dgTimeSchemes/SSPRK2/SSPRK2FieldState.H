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
    Foam::SSPRK2FieldState

Description
    Storage container for one DG field advanced by the two-stage SSPRK2
    time scheme.

SourceFiles
    SSPRK2FieldStateI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_SSPRK2FieldState_H
#define Foam_SSPRK2FieldState_H

#include "dgField.H"
#include "dgGeomMesh.H"
#include "IOobject.H"
#include "List.H"
#include "PtrList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SSPRK2FieldState Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::SSPRK2FieldState
 * \brief Per-field state and residual storage for SSPRK2.
 *
 * \tparam Type DG field value type, typically \c scalar or \c vector
 *
 * \details
 * The owning time-scheme applies
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
 * To support this update, the state stores the current field \c U, the
 * beginning-of-step snapshot \c UOld, one persistent intermediate field
 * \c UStage1, and two residual buffers. The object is intentionally
 * lightweight in behavior and acts mainly as a structured container used by
 * \c dgSSPRK2TimeScheme. All stage algebra is performed in the owning scheme;
 * this class only manages field snapshots and residual buffers with
 * mesh-consistent sizing.
 */
template<class Type>
class SSPRK2FieldState
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Number of SSPRK2 residual stages. */
        static constexpr label nStages_ = 2;

        /** Number of intermediate stage fields that must be stored. */
        static constexpr label nStoredStages_ = 1;

        /** DG mesh used to size per-cell residual and stage storage. */
        const dgGeomMesh& mesh_;

        /** Active conservative DG field currently being advanced. */
        dgField<Type>& U_;

        /** Snapshot of \f$\mathbf{U}^{n}\f$ taken at the start of the time step. */
        dgField<Type> UOld_;

        /** Stored stage field, currently representing \f$\mathbf{U}^{(1)}\f$. */
        PtrList<dgField<Type>> UStages_;

        /** Residual buffers \f$\mathbf{R}_s\f$ indexed as [stage][cell][dof]. */
        List<List<List<Type>>> RStages_;


    //====================================================================//
    //                       Internal helper functions
    //====================================================================//

        /** \brief Allocate intermediate stage fields with the same layout as \c U_. */
        void allocateStageFields();

        /** \brief Allocate residual buffers using the local number of DoFs per cell. */
        void allocateResiduals();

        /**
         * \brief Copy one DG field into another.
         *
         * \param dst Destination field
         * \param src Source field
         */
        void copyField(dgField<Type>& dst, const dgField<Type>& src);

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** Disallow default construction. */
        SSPRK2FieldState() = delete;

        /**
         * \brief Construct state storage for one registered field.
         *
         * \param field Active DG field to be managed
         * \param mesh DG mesh shared by the solver
         */
        SSPRK2FieldState
        (
            dgField<Type>& field,
            const dgGeomMesh& mesh
        );

        /** Default destructor. */
        ~SSPRK2FieldState() = default;


    //====================================================================//
    //                              Access
    //====================================================================//

        /** \brief Return the active field being advanced. */
        inline dgField<Type>& U()
        {
            return U_;
        }

        /** \brief Return the active field being advanced. */
        inline const dgField<Type>& U() const
        {
            return U_;
        }

        /** \brief Return the beginning-of-step field snapshot \f$\mathbf{U}^{n}\f$. */
        inline dgField<Type>& UOld()
        {
            return UOld_;
        }

        /** \brief Return the beginning-of-step field snapshot \f$\mathbf{U}^{n}\f$. */
        inline const dgField<Type>& UOld() const
        {
            return UOld_;
        }

        /** \brief Return the stored first stage field \f$\mathbf{U}^{(1)}\f$. */
        inline dgField<Type>& UStage1()
        {
            return UStages_[0];
        }

        /** \brief Return the stored first stage field \f$\mathbf{U}^{(1)}\f$. */
        inline const dgField<Type>& UStage1() const
        {
            return UStages_[0];
        }

        /**
         * \brief Return one stored intermediate stage field by index.
         *
         * \param stageFieldI Zero-based stage-field index
         */
        inline dgField<Type>& UStage(const label stageFieldI)
        {
            return UStages_[stageFieldI];
        }

        /**
         * \brief Return one stored intermediate stage field by index.
         *
         * \param stageFieldI Zero-based stage-field index
         */
        inline const dgField<Type>& UStage(const label stageFieldI) const
        {
            return UStages_[stageFieldI];
        }

        /** \brief Return the number of SSPRK2 residual stages. */
        inline label nStages() const
        {
            return nStages_;
        }

        /** \brief Return the number of stored intermediate stage fields. */
        inline label nStoredStages() const
        {
            return nStoredStages_;
        }

        /** \brief Return all residual buffers. */
        inline List<List<List<Type>>>& RStages()
        {
            return RStages_;
        }

        /** \brief Return all residual buffers. */
        inline const List<List<List<Type>>>& RStages() const
        {
            return RStages_;
        }

        /**
         * \brief Return the residual buffer of one stage.
         *
         * \param stageI Zero-based SSPRK2 stage index
         */
        inline List<List<Type>>& R(const label stageI)
        {
            return RStages_[stageI];
        }

        /**
         * \brief Return the residual buffer of one stage.
         *
         * \param stageI Zero-based SSPRK2 stage index
         */
        inline const List<List<Type>>& R(const label stageI) const
        {
            return RStages_[stageI];
        }

        /**
         * \brief Return the residual vector for one cell and stage.
         *
         * \param stageI Zero-based SSPRK2 stage index
         * \param cellI Local DG cell index
         */
        inline List<Type>& R(const label stageI, const label cellI)
        {
            return RStages_[stageI][cellI];
        }

        /**
         * \brief Return the residual vector for one cell and stage.
         *
         * \param stageI Zero-based SSPRK2 stage index
         * \param cellI Local DG cell index
         */
        inline const List<Type>& R(const label stageI, const label cellI) const
        {
            return RStages_[stageI][cellI];
        }


    //====================================================================//
    //                              State API
    //====================================================================//

        /** \brief Save the active field into \c UOld_. */
        void saveOldField();

        /** \brief Restore the active field from \c UOld_. */
        void restoreOldField();

        /**
         * \brief Save the active field into one stored intermediate stage field.
         *
         * \param stageFieldI Zero-based stage-field index
         */
        void saveStageField(const label stageFieldI);

        /**
         * \brief Restore the active field from one stored intermediate stage field.
         *
         * \param stageFieldI Zero-based stage-field index
         */
        void restoreStageField(const label stageFieldI);

        /** \brief Set every residual entry in every stage buffer to zero. */
        void clearResiduals();

        /**
         * \brief Set every residual entry in one stage buffer to zero.
         *
         * \param stageI Zero-based SSPRK2 stage index
         */
        void clearResidual(const label stageI);

        /** \brief Restore a clean state by copying back \c UOld_ and clearing residuals. */
        void reset();
};

} // End namespace Foam

#include "SSPRK2FieldStateI.H"

#endif

// ************************************************************************* //
