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
    Foam::SSPRK3FieldState

Description
    Storage container for one DG field advanced by the three-stage SSPRK3
    or TVDRK3 time scheme.

SourceFiles
    SSPRK3FieldStateI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_SSPRK3FieldState_H
#define Foam_SSPRK3FieldState_H

#include "dgField.H"
#include "dgGeomMesh.H"
#include "IOobject.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SSPRK3FieldState Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::SSPRK3FieldState
 * \brief Per-field state and residual storage for SSPRK3 / TVDRK3.
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
 * \mathbf{U}^{(2)} &=
 *     \frac{3}{4}\mathbf{U}^{n}
 *   + \frac{1}{4}
 *     \left(
 *         \mathbf{U}^{(1)}
 *       + \Delta t\,\mathbf{M}^{-1}\mathbf{R}_2
 *     \right), \\
 * \mathbf{U}^{n+1} &=
 *     \frac{1}{3}\mathbf{U}^{n}
 *   + \frac{2}{3}
 *     \left(
 *         \mathbf{U}^{(2)}
 *       + \Delta t\,\mathbf{M}^{-1}\mathbf{R}_3
 *     \right).
 * \end{aligned}
 * \f]
 *
 * The object stores the active field \c U, the beginning-of-step snapshot
 * \c UOld, and three residual buffers. Because every stage reuses only the
 * current solution and the saved beginning-of-step field, no persistent
 * intermediate stage field is needed.
 */
template<class Type>
class SSPRK3FieldState
{
private:

        /** Number of SSPRK3 residual stages. */
        static constexpr label nStages_ = 3;

        /** DG mesh used to size per-cell residual storage. */
        const dgGeomMesh& mesh_;

        /** Active conservative DG field currently being advanced. */
        dgField<Type>& U_;

        /** Snapshot of \f$\mathbf{U}^{n}\f$ taken at the start of the time step. */
        dgField<Type> UOld_;

        /** Residual buffers \f$\mathbf{R}_s\f$ indexed as [stage][cell][dof]. */
        List<List<List<Type>>> RStages_;

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

        /** Disallow default construction. */
        SSPRK3FieldState() = delete;

        /**
         * \brief Construct state storage for one registered field.
         *
         * \param field Active DG field to be managed
         * \param mesh DG mesh shared by the solver
         */
        SSPRK3FieldState
        (
            dgField<Type>& field,
            const dgGeomMesh& mesh
        );

        /** Default destructor. */
        ~SSPRK3FieldState() = default;

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

        /** \brief Return the number of SSPRK3 residual stages. */
        inline label nStages() const
        {
            return nStages_;
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
         * \param stageI Zero-based SSPRK3 stage index
         */
        inline List<List<Type>>& R(const label stageI)
        {
            return RStages_[stageI];
        }

        /**
         * \brief Return the residual buffer of one stage.
         *
         * \param stageI Zero-based SSPRK3 stage index
         */
        inline const List<List<Type>>& R(const label stageI) const
        {
            return RStages_[stageI];
        }

        /**
         * \brief Return the residual vector for one cell and stage.
         *
         * \param stageI Zero-based SSPRK3 stage index
         * \param cellI Local DG cell index
         */
        inline List<Type>& R(const label stageI, const label cellI)
        {
            return RStages_[stageI][cellI];
        }

        /**
         * \brief Return the residual vector for one cell and stage.
         *
         * \param stageI Zero-based SSPRK3 stage index
         * \param cellI Local DG cell index
         */
        inline const List<Type>& R(const label stageI, const label cellI) const
        {
            return RStages_[stageI][cellI];
        }

        /** \brief Save the active field into \c UOld_. */
        void saveOldField();

        /** \brief Restore the active field from \c UOld_. */
        void restoreOldField();

        /** \brief Set every residual entry in every stage buffer to zero. */
        void clearResiduals();

        /**
         * \brief Set every residual entry in one stage buffer to zero.
         *
         * \param stageI Zero-based SSPRK3 stage index
         */
        void clearResidual(const label stageI);

        /** \brief Restore a clean state by copying back \c UOld_ and clearing residuals. */
        void reset();
};

} // End namespace Foam

#include "SSPRK3FieldStateI.H"

#endif

// ************************************************************************* //
