/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Store degrees of freedom (Dof) of given Type for a DG cell.

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

template<class Type>
class cellDof
{
private:

    //- Cell ID in the mesh
    label cellID_;

    //- Number of degrees of freedom stored in this cell
    label nDof_;

    //- List of DoF values
    List<Type> dof_;

public:

    // Constructors

        //- Default constructor (required for List<T> containers)
        cellDof();

        //- Copy constructor
        cellDof(const cellDof<Type>& other);

        //- Construct from cell ID and an input DoF list
        cellDof
        (
            const label cellID,
            const label nDof,
            const UList<Type>& inputDof
        );

    //- Assignment operator
    cellDof<Type>& operator=(const cellDof<Type>& other);

    // Access

        //- Return non-const access to DoF list
        inline List<Type>& dof()
        {
            return dof_;
        }

        //- Return const access to DoF list
        inline const List<Type>& dof() const
        {
            return dof_;
        }

        //- Return number of DoF in this cell
        inline label nDof() const
        {
            return nDof_;
        }

        //- Return cell ID
        inline label cellID() const
        {
            return cellID_;
        }

        //- Non-const element access
        inline Type& operator[](const label i)
        {
            return dof_[i];
        }

        //- Const element access
        inline const Type& operator[](const label i) const
        {
            return dof_[i];
        }
};

} // End namespace Foam

#include "cellDofI.H"

#endif

// ************************************************************************* //

