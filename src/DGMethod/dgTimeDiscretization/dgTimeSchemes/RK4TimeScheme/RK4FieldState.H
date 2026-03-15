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
    Foam::RK4FieldState

Description
    Per-field storage used by the classical explicit RK4 time scheme.

SourceFiles
    RK4FieldStateI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_RK4FieldState_H
#define Foam_RK4FieldState_H

#include "dgField.H"
#include "dgGeomMesh.H"
#include "IOobject.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class RK4FieldState Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::RK4FieldState
 * \brief Per-field storage used by the classical four-stage Runge-Kutta scheme.
 *
 * \details
 * For the classical RK4 update,
 *
 * \f[
 *     \mathbf{U}^{n+1}
 *     =
 *     \mathbf{U}^{n}
 *     +
 *     \frac{\Delta t}{6}
 *     \left(
 *         \mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4
 *     \right),
 * \f]
 *
 * the scheme stores the beginning-of-step state \f$\mathbf{U}^n\f$ together
 * with four residual buffers corresponding to the stage vectors
 * \f$\mathbf{k}_1,\dots,\mathbf{k}_4\f$.
 */
template<class Type>
class RK4FieldState
{
private:

        /** Number of stages in classical RK4. */
        static constexpr label nStages_ = 4;

        /** Shared DG mesh used to size stage residual buffers. */
        const dgGeomMesh& mesh_;

        /** Active DG field currently being advanced in time. */
        dgField<Type>& U_;

        /** Snapshot of the beginning-of-step field \f$\mathbf{U}^n\f$. */
        dgField<Type> UOld_;

        /** Stage residual storage indexed as [stage][cell][dof]. */
        List<List<List<Type>>> RStages_;

        /** \brief Allocate residual buffers according to the local cell DoF counts. */
        void allocateResiduals();

        /**
         * \brief Copy one DG field into another internal storage field.
         *
         * \param dst Destination field
         * \param src Source field
         */
        void copyField(dgField<Type>& dst, const dgField<Type>& src);

public:

        /** Disallow default construction. */
        RK4FieldState() = delete;

        /**
         * \brief Construct RK4 state for one DG field.
         *
         * \param field Active DG field
         * \param mesh DG mesh associated with the field
         */
        RK4FieldState
        (
            dgField<Type>& field,
            const dgGeomMesh& mesh
        );

        /** Default destructor. */
        ~RK4FieldState() = default;

        /** \brief Return writable access to the active DG field. */
        inline dgField<Type>& U()
        {
            return U_;
        }

        /** \brief Return read-only access to the active DG field. */
        inline const dgField<Type>& U() const
        {
            return U_;
        }

        /** \brief Return writable access to the saved beginning-of-step field. */
        inline dgField<Type>& UOld()
        {
            return UOld_;
        }

        /** \brief Return read-only access to the saved beginning-of-step field. */
        inline const dgField<Type>& UOld() const
        {
            return UOld_;
        }

        /** \brief Return the number of RK4 stages, namely 4. */
        inline label nStages() const
        {
            return nStages_;
        }

        /** \brief Return writable access to the complete stage-residual storage. */
        inline List<List<List<Type>>>& RStages()
        {
            return RStages_;
        }

        /** \brief Return read-only access to the complete stage-residual storage. */
        inline const List<List<List<Type>>>& RStages() const
        {
            return RStages_;
        }

        /**
         * \brief Return the residual buffer of one RK4 stage.
         *
         * \param stageI Zero-based stage index
         * \return Writable residual storage indexed by [cell][dof]
         */
        inline List<List<Type>>& R(const label stageI)
        {
            return RStages_[stageI];
        }

        /**
         * \brief Return the residual buffer of one RK4 stage.
         *
         * \param stageI Zero-based stage index
         * \return Read-only residual storage indexed by [cell][dof]
         */
        inline const List<List<Type>>& R(const label stageI) const
        {
            return RStages_[stageI];
        }

        /**
         * \brief Return the residual vector of one RK4 stage and one cell.
         *
         * \param stageI Zero-based stage index
         * \param cellI Local DG cell index
         * \return Writable residual vector indexed by DoF
         */
        inline List<Type>& R(const label stageI, const label cellI)
        {
            return RStages_[stageI][cellI];
        }

        /**
         * \brief Return the residual vector of one RK4 stage and one cell.
         *
         * \param stageI Zero-based stage index
         * \param cellI Local DG cell index
         * \return Read-only residual vector indexed by DoF
         */
        inline const List<Type>& R(const label stageI, const label cellI) const
        {
            return RStages_[stageI][cellI];
        }

        /** \brief Return the number of cells covered by this field state. */
        inline label nCells() const
        {
            return U_.nCells();
        }

        /** \brief Save the current field as the beginning-of-step state. */
        void saveOldField();

        /** \brief Restore the active field from the saved beginning-of-step state. */
        void restoreOldField();

        /** \brief Set all stage residual buffers to zero. */
        void clearResiduals();

        /**
         * \brief Set one stage residual buffer to zero.
         *
         * \param stageI Zero-based stage index
         */
        void clearResidual(const label stageI);

        /** \brief Restore the old field and clear all stage buffers. */
        void reset();
};

} // End namespace Foam

#include "RK4FieldStateI.H"

#endif

// ************************************************************************* //
