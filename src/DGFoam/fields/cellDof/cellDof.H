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

    //- Cell ID
    label cellID_;

    //- Number of degrees of freedom
    label nDof_;

    //- List of Dof values of Type
    List<Type> dof_;

public:

    //- Default constructor (required for List<T>)
    cellDof();

    //- Copy constructor
    cellDof(const cellDof<Type>& other);

    //- Construct from cell ID, number of DoF, and input Dof list
    cellDof
    (
        const label cellID,
        const label nDof,
        const UList<Type>& inputDof
    );

    //- Assignment operator
    cellDof<Type>& operator=(const cellDof<Type>& other);

    // Accessors

    inline List<Type>& dof() { return dof_; }

    inline const List<Type>& dof() const { return dof_; }

    inline label nDof() const { return nDof_; }

    inline label cellID() const { return cellID_; }

    inline Type& operator[](const label i) { return dof_[i]; }

    inline const Type& operator[](const label i) const { return dof_[i]; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
