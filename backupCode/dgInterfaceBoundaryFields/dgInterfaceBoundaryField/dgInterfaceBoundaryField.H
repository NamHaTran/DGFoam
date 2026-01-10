/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DGFoam: Discontinuous Galerkin CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | GPU-friendly CFD solver framework
     \\/     M anipulation  |
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgInterfaceBoundaryField

Description
    Base class for runtime-selectable DG boundary conditions.
    When creating a new boundary condition, please remember to override
    checkPatchType() method.

    For dgInterfaceBoundaryField, the runtime selection is based on
    fvPatch::type(), not on a dictionary entry.
    Therefore, all derived dgInterfaceBoundaryField classes must have
    TypeName equal to the associated fvPatch type (e.g. "processor",
    "cyclic", "cyclicAMI").
    This is fundamentally different from dgGeneralBoundaryField, where
    the selection key comes from the "type" entry in the boundary
    condition dictionary.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgInterfaceBoundaryField_H
#define Foam_dgInterfaceBoundaryField_H

#include "word.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dgGeomMesh.H"
#include "fvPatch.H"
#include "dgField.H"

namespace Foam
{

// Forward declaration
template<class Type> class dgField;

/*---------------------------------------------------------------------------*\
                    Class dgInterfaceBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgInterfaceBoundaryField
{
protected:

    //- Reference to fvPatch (geometry of patch)
    const fvPatch& patch_;

    //- Reference to DG geometric mesh (for DG-level cell/face connectivity)
    const dgGeomMesh& dgMesh_;

public:

    //- Runtime type information
    TypeName("dgInterfaceBoundaryField");

    //- Runtime selection table (dictionary constructor)
    declareRunTimeSelectionTable
    (
        autoPtr,
        dgInterfaceBoundaryField,
        fvPatch,
        (const fvPatch& patch, const dgGeomMesh& dgMesh),
        (patch, dgMesh)
    );

    // Constructors

        //- Construct from fvPatch and dgGeomMesh (default dictionary)
        dgInterfaceBoundaryField
        (
            const fvPatch& patch,
            const dgGeomMesh& dgMesh
        );

    //- Destructor
    virtual ~dgInterfaceBoundaryField() = default;

    // Factory Method

        static autoPtr<dgInterfaceBoundaryField<Type>> New
        (
            const fvPatch& patch,
            const dgGeomMesh& dgMesh
        );

    // Core Interface

        virtual void update
        (
            dgField<Type>& field
        ) const
        { }

    // Access Functions

        inline const fvPatch& patch() const { return patch_; }

        inline const dgGeomMesh& dgMesh() const { return dgMesh_; }

        inline const word& name() const { return patch_.name(); }

        inline word BCtype() const
        {
            return patch_.type();
        }

    // Check type
        
        //- Verify patch type compatibility
        virtual void checkPatchType() const {}
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgInterfaceBoundaryFieldI.H"

#endif
// ************************************************************************* //
