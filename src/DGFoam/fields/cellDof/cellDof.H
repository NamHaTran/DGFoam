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
    Container for degrees of freedom (DoF) of a DG cell.

    Stores all polynomial degrees of freedom associated with a single
    DG cell for a given field type.

SourceFiles
    cellDof.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellDof_H
#define Foam_cellDof_H

#include "label.H"
#include "List.H"

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
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Cell index in the mesh */
        label cellID_;

        /** Number of degrees of freedom stored in this cell */
        label nDof_;

        /** List of degree-of-freedom values */
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
         * \brief Copy constructor
         */
        cellDof(const cellDof<Type>& other);

        /**
         * \brief Construct from cell ID and DoF values
         *
         * Copies the provided DoF values into internal storage.
         */
        cellDof
        (
            const label cellID,
            const label nDof,
            const UList<Type>& inputDof
        );

    //====================================================================//
    //                         Assignment
    //====================================================================//

        /**
         * \brief Copy assignment operator
         */
        cellDof<Type>& operator=(const cellDof<Type>& other);

    //====================================================================//
    //                           Access
    //====================================================================//

        /** Return non-const access to DoF list */
        inline List<Type>& dof()
        {
            return dof_;
        }

        /** Return const access to DoF list */
        inline const List<Type>& dof() const
        {
            return dof_;
        }

        /** Return number of DoFs in this cell */
        inline label nDof() const
        {
            return nDof_;
        }

        /** Return cell index */
        inline label cellID() const
        {
            return cellID_;
        }

        /** Non-const access to i-th DoF */
        inline Type& operator[](const label i)
        {
            return dof_[i];
        }

        /** Const access to i-th DoF */
        inline const Type& operator[](const label i) const
        {
            return dof_[i];
        }
};

} // End namespace Foam

#include "cellDofI.H"

#endif

// ************************************************************************* //
