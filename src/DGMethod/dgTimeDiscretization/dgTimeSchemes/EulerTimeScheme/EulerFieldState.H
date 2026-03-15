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
    Foam::EulerFieldState

Description
    Per-field storage used by the explicit Euler time scheme.

    EulerFieldState stores exactly the data needed by one-step explicit Euler:
    - a reference to the active field \f$\mathbf{U}\f$
    - a saved copy of the beginning-of-step state \f$\mathbf{U}^n\f$
    - one residual buffer \f$\mathbf{R}^n\f$ organized by [cell][dof]

SourceFiles
    EulerFieldStateI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_EulerFieldState_H
#define Foam_EulerFieldState_H

#include "dgField.H"
#include "dgGeomMesh.H"
#include "IOobject.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class EulerFieldState Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::EulerFieldState
 * \brief Per-field storage used by explicit Euler time integration.
 *
 * \details
 * For the explicit Euler update
 *
 * \f[
 *     \mathbf{U}^{n+1}
 *     =
 *     \mathbf{U}^{n}
 *     +
 *     \Delta t\,\mathbf{M}^{-1}\mathbf{R}^{n},
 * \f]
 *
 * the scheme needs only one saved field \f$\mathbf{U}^n\f$ and one residual
 * buffer \f$\mathbf{R}^{n}\f$ per registered DG field.
 */
template<class Type>
class EulerFieldState
{
private:

        /** Shared DG mesh used to size residual buffers from the local DoF count. */
        const dgGeomMesh& mesh_;

        /** Active DG field currently being advanced in time. */
        dgField<Type>& U_;

        /** Snapshot of the beginning-of-step field \f$\mathbf{U}^n\f$. */
        dgField<Type> UOld_;

        /** Residual storage indexed as [cell][dof]. */
        List<List<Type>> R_;

        /**
         * \brief Allocate the residual storage according to the cell-local number of DoFs.
         */
        void allocateResidual();

        /**
         * \brief Copy one DG field into another internal storage field.
         *
         * \param dst Destination field
         * \param src Source field
         */
        void copyField(dgField<Type>& dst, const dgField<Type>& src);

public:

        /** Disallow default construction. */
        EulerFieldState() = delete;

        /**
         * \brief Construct Euler state for one DG field.
         *
         * The saved field uses \c NO_READ/\c NO_WRITE because it is purely an
         * internal time-integration buffer.
         *
         * \param field Active DG field
         * \param mesh DG mesh associated with the field
         */
        EulerFieldState
        (
            dgField<Type>& field,
            const dgGeomMesh& mesh
        );

        /** Default destructor. */
        ~EulerFieldState() = default;

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

        /** \brief Return writable access to the whole residual buffer. */
        inline List<List<Type>>& R()
        {
            return R_;
        }

        /** \brief Return read-only access to the whole residual buffer. */
        inline const List<List<Type>>& R() const
        {
            return R_;
        }

        /**
         * \brief Return the residual vector of one cell.
         *
         * \param cellI Local DG cell index
         * \return Writable residual vector indexed by DoF
         */
        inline List<Type>& R(const label cellI)
        {
            return R_[cellI];
        }

        /**
         * \brief Return the residual vector of one cell.
         *
         * \param cellI Local DG cell index
         * \return Read-only residual vector indexed by DoF
         */
        inline const List<Type>& R(const label cellI) const
        {
            return R_[cellI];
        }

        /** \brief Return the number of cells covered by this field state. */
        inline label nCells() const
        {
            return R_.size();
        }

        /** \brief Save the current field as the beginning-of-step state. */
        void saveOldField();

        /** \brief Restore the active field from the saved beginning-of-step state. */
        void restoreOldField();

        /** \brief Set the complete residual buffer to zero. */
        void clearResidual();

        /** \brief Restore the old field and clear the residual buffer. */
        void reset();
};

} // End namespace Foam

#include "EulerFieldStateI.H"

#endif

// ************************************************************************* //
