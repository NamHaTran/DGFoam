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
    Foam::dgField

Description
    DG field class storing DOF and GaussField data for the whole mesh.

SourceFiles
    dgBasisField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgField_H
#define Foam_dgField_H

#include "regIOobject.H"
#include "dofField.H"
#include "GaussField.H"
#include "dgGeomMesh.H"
#include "dgGeneralBoundaryManager.H"
#include "tensor.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class dgField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgField
:
    public regIOobject
{
private:

    //- Reference to DG geometric mesh
    const dgGeomMesh& mesh_;

    //- True if this field owns DOF data
    bool hasDof_;

    //- Pointer to DOF field (null for primary fields)
    autoPtr<dofField<Type>> dofPtr_;

    //- GaussField for every cell in the mesh
    List<GaussField<Type>> gaussFields_;

public:

    //- Runtime type information
    TypeName("dgField");

    // * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * * //

    dgField
    (
        const IOobject& io,
        const dgGeomMesh& mesh,
        bool hasDof = true
    );

    //- No copy
    dgField(const dgField&) = delete;
    dgField& operator=(const dgField&) = delete;

    // * * * * * * * * * * * Member Functions * * * * * * * * * * * //

    //- Return true if this field has DOF data
    inline bool hasDof() const
    {
        return dofPtr_ != nullptr;
    }

    //- Return non-const access to DOF field
    inline dofField<Type>& dof()
    {
        if (!dofPtr_)
        {
            FatalErrorInFunction
                << "Attempt to access non-existing DOF field"
                << abort(FatalError);
        }
        return *dofPtr_;
    }

    //- Return const access to DOF field
    inline const dofField<Type>& dof() const
    {
        if (!dofPtr_)
        {
            FatalErrorInFunction
                << "Attempt to access non-existing DOF field"
                << abort(FatalError);
        }
        return *dofPtr_;
    }

    //- Return non-const access to GaussFields
    inline List<GaussField<Type>>& gaussFields()
    {
        return gaussFields_;
    }

    //- Return const access to GaussFields
    inline const List<GaussField<Type>>& gaussFields() const
    {
        return gaussFields_;
    }

    inline label nCells() const
    {
        return gaussFields_.size();
    } 

    //- Write data (minimal)
    bool writeData(Ostream& os) const;

    // Assignment ------------------------------------------------------------
    
    dgField<Type>& operator=(const Type& value);
};

} // End namespace Foam

#include "dgFieldI.H"

#endif

// ************************************************************************* //
