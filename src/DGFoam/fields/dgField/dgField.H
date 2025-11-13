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
#include "autoPtr.H"
#include "dofField.H"
#include "GaussField.H"
#include "dgGeomMesh.H"
#include "dgGeneralBoundaryManager.H"
#include "tensor.H"
#include "dofField.H"

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

    //- Reference to the DG geometric mesh
    const dgGeomMesh& mesh_;

    //- Has DOF field
    bool hasDof_;

    //- Pointer to DOF field (only valid for conservative fields)
    // autoPtr<dofField<Type>> dofPtr_;
    dofField<Type>* dofPtr_ = nullptr;

    //- List of Gauss fields for all cells in the mesh
    List<GaussField<Type>> gaussFields_;

    //- Pointer to boundary manager (optional)
    dgGeneralBoundaryManager<Type>* bcManagerPtr_;

public:

    //- Runtime type information
    TypeName("dgField");

    // Constructors
    //-------------------------------------------------------------------------

    //- Construct from IOobject and mesh
    //  If hasDof == true, the field will read DOF from file and interpolate Gauss values.
    //  Otherwise, the field will create empty GaussField structures (for primary fields).
    dgField
    (
        const IOobject& io,
        const dgGeomMesh& mesh,
        bool hasDof = true
    );

    //- Disable copy constructor and assignment operator
    dgField(const dgField&) = delete;
    dgField& operator=(const dgField&) = delete;

    //- Destructor
    ~dgField();


    // Member Functions
    //-------------------------------------------------------------------------

    //- Return true if DOF field is available
    inline bool hasDof() const
    {
        return dofPtr_ != nullptr;
    }

    //- Access to DOF field (non-const)
    inline dofField<Type>& dof()
    {
        if (!dofPtr_)
        {
            FatalErrorInFunction
                << "Attempting to access DOF field in dgField<"
                << pTraits<Type>::typeName
                << "> but no DOF data has been created for this field."
                << abort(FatalError);
        }

        return *dofPtr_;
    }

    //- Access to DOF field (const)
    inline const dofField<Type>& dof() const
    {
        if (!dofPtr_)
        {
            FatalErrorInFunction
                << "Attempting to access DOF field in dgField<"
                << pTraits<Type>::typeName
                << "> but no DOF data has been created for this field."
                << abort(FatalError);
        }

        return *dofPtr_;
    }

    //- Return list of Gauss fields (non-const)
    inline List<GaussField<Type>>& gaussFields()
    {
        return gaussFields_;
    }

    //- Return list of Gauss fields (const)
    inline const List<GaussField<Type>>& gaussFields() const
    {
        return gaussFields_;
    }

    //- Return Gauss field of a specific cell
    inline GaussField<Type>& cellGaussField(const label cellI)
    {
        return gaussFields_[cellI];
    }

    //- Return number of cells
    inline label nCells() const
    {
        return gaussFields_.size();
    }

    //- Set boundary manager pointer
    inline void setBCManager(dgGeneralBoundaryManager<Type>* ptr)
    {
        bcManagerPtr_ = ptr;
    }

    //- Return boundary manager pointer
    inline dgGeneralBoundaryManager<Type>* bcManager() const
    {
        return bcManagerPtr_;
    }

    //- Placeholder for boundary correction (to be implemented later)
    void correctBoundary();

    //- Write data
    bool writeData(Ostream& os) const;

    // I/O
    //-------------------------------------------------------------------------
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
