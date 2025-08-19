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
#include "dgGeomMesh.H"
#include "dgGeomCell.H"

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

    //- Reference to DG geometry mesh
    dgGeomMesh& dgMesh_;

    //- Reference to geometric cell
    dgGeomCell& cell_;

    //- Number of degrees of freedom
    label nDof_;

    //- List of Dof values of Type
    List<Type> dof_;

public:

    //- Construct from cell ID, mesh, and input Dof list
    cellDof
    (
        const label cellID,
        dgGeomMesh& dgMesh,
        const UList<Type>& inputDof
    );

    //- Return modifiable reference to Dof list
    inline List<Type>& dof()
    {
        return dof_;
    }

    //- Return const reference to Dof list
    inline const List<Type>& dof() const
    {
        return dof_;
    }

    //- Return number of Dof
    inline label nDof() const
    {
        return nDof_;
    }

    //- Return cell ID
    inline label cellID() const
    {
        return cellID_;
    }

    //- Return reference to DG mesh
    inline dgGeomMesh& mesh() const
    {
        return dgMesh_;
    }

    //- Return reference to DG geometric cell
    inline const dgGeomCell& cell() const
    {
        return cell_;
    }

    //- Return non-const reference to Dof at given index
    inline Type& operator[](const label i)
    {
        return dof_[i];
    }

    //- Return const reference to Dof at given index
    inline const Type& operator[](const label i) const
    {
        return dof_[i];
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
