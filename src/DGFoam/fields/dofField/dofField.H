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
    Foam::dofField<Type>

Description
    Store degrees of freedom (Dof) of given Type for a DG cell.

SourceFiles
    dofField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dofField_H
#define Foam_dofField_H

#include "volFields.H"
#include "cellDof.H"
#include "dgGeomMesh.H"
#include "IOobject.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dofField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dofField
{
private:

    //- Reference to DG geometric mesh
    const dgGeomMesh& dgMesh_;

    //- List of cellDoF objects (local DoF storage)
    List<cellDof<Type>> cellDofs_;

    //- Number of cells
    label nCells_;

    //- List of OpenFOAM volFields, one field per DoF index
    List<autoPtr<GeometricField<Type, fvPatchField, volMesh>>> foamFields_;

    //- Number of DoFs for each cell
    List<label> dofPerCell_;

    //- Maximum DoF across all cells
    label maxDoF_;

public:

    // Constructors ------------------------------------------------------------

    //- Construct from field base name, mesh and read/write options
    dofField
    (
        const word& name,
        const fileName& instance,
        IOobject::readOption rOpt,
        IOobject::writeOption wOpt,
        const dgGeomMesh& mesh
    );

    // Access -----------------------------------------------------------------

    inline const List<cellDof<Type>>& cellDofList() const
    {
        return cellDofs_;
    }

    inline List<cellDof<Type>>& cellDofList()
    {
        return cellDofs_;
    }

    inline cellDof<Type>& operator[](const label i)
    {
        return cellDofs_[i];
    }

    inline const cellDof<Type>& operator[](const label i) const
    {
        return cellDofs_[i];
    }

    // Member Functions --------------------------------------------------------

    //- Update OpenFOAM volFields from DG DoF values
    void updateFoamFields();

    //- Write volFields to disk
    void write() const;
};

} // End namespace Foam

#include "dofFieldI.H"

#endif

// ************************************************************************* //
