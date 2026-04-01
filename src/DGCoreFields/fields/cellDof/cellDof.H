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

Class
    Foam::cellDof<Type>

Description
    Degree-of-freedom container for a DG cell.

    Stores all polynomial degrees of freedom associated with a single
    DG cell for a given field type.

SourceFiles
    cellDofI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellDof_H
#define Foam_cellDof_H

#include "label.H"
#include "List.H"
#include "tmp.H"
#include "refCount.H"

// Forward declarations of tensor types
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellDof Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::cellDof
 * \brief Degree-of-freedom container for a DG cell
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * cellDof stores all degrees of freedom (DoF) associated with a single
 * DG cell for a given field type.
 *
 * The DoFs typically correspond to coefficients of polynomial basis
 * functions defined on the reference cell and mapped to physical space.
 *
 * This class is intentionally lightweight and contains no geometric
 * information; it is purely a data container used by higher-level DG
 * field and solver classes.
 */
template<class Type>
class cellDof
:
    public refCount
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** \brief Index of the mesh cell represented by this DoF container. */
        label cellID_;

        /** \brief Number of degrees of freedom stored for the cell. */
        label nDof_;

        /** \brief Polynomial degree-of-freedom values ordered by local basis index. */
        List<Type> dof_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Default constructor
         *
         * Required for List<cellDof<Type>> containers.
         */
        cellDof();

        /**
         * \brief Copy constructor.
         *
         * \param other Source cell DoF container to copy from.
         */
        cellDof(const cellDof<Type>& other);

        /**
         * \brief Construct from a cell index and an external DoF list.
         *
         * Copies up to \p nDof values from \p inputDof into internal storage.
         * Remaining entries, if any, are initialised to \c pTraits<Type>::zero.
         *
         * \param cellID Mesh cell index associated with the DoF set.
         * \param nDof Number of DoFs to allocate for the cell.
         * \param inputDof External DoF values copied into this object.
         */
        cellDof
        (
            const label cellID,
            const label nDof,
            const UList<Type>& inputDof
        );

        /**
         * \brief Construct by transferring storage from a temporary cellDof.
         *
         * \param t Temporary cell DoF container whose storage is consumed.
         */
        cellDof
        (
            const tmp<cellDof<Type>>& t
        );

    //====================================================================//
    //                         Assignment
    //====================================================================//

        /**
         * \brief Copy assignment operator.
         *
         * \param other Source cell DoF container to copy from.
         * \return Reference to this object.
         */
        cellDof<Type>& operator=(const cellDof<Type>& other);

        /**
         * \brief Assign by transferring storage from a temporary cellDof.
         *
         * \param t Temporary cell DoF container whose storage is consumed.
         * \return Reference to this object.
         */
        cellDof<Type>& operator=
        (
            const tmp<cellDof<Type>>& t
        );

    //====================================================================//
    //                           Access
    //====================================================================//

        /** \brief Return writable access to the complete DoF storage. */
        inline List<Type>& dof()
        {
            return dof_;
        }

        /** \brief Return read-only access to the complete DoF storage. */
        inline const List<Type>& dof() const
        {
            return dof_;
        }

        /** \brief Return the number of degrees of freedom stored in the cell. */
        inline label nDof() const
        {
            return nDof_;
        }

        /** \brief Return the mesh cell index associated with this DoF set. */
        inline label cellID() const
        {
            return cellID_;
        }

        /**
         * \brief Return writable access to one local degree of freedom.
         *
         * \param i Local DoF index inside the cell.
         * \return Reference to the requested DoF value.
         */
        inline Type& operator[](const label i)
        {
            return dof_[i];
        }

        /**
         * \brief Return read-only access to one local degree of freedom.
         *
         * \param i Local DoF index inside the cell.
         * \return Const reference to the requested DoF value.
         */
        inline const Type& operator[](const label i) const
        {
            return dof_[i];
        }
};

} // End namespace Foam

#include "cellDofI.H"

#endif

// ************************************************************************* //
