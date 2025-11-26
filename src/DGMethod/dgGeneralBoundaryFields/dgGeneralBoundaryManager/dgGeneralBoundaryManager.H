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
    Foam::dgGeneralBoundaryManager

Description
    Class for DG boundary conditions management.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeneralBoundaryManager_H
#define Foam_dgGeneralBoundaryManager_H

#include "dgGeomMesh.H"
#include "fvPatch.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "dgGeneralBoundaryField.H"
#include "dgField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgGeneralBoundaryManager Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgGeneralBoundaryManager
{
private:

    //- Reference to dgGeomMesh (used for fvPatch access and field lookup)
    const dgGeomMesh& mesh_;

    //- List of boundary condition objects (one per fvPatch)
    PtrList<dgGeneralBoundaryField<Type>> bConditions_;

    //- Internal field value (currently only supports 'uniform <value>')
    Type internalValue_;

public:

    // Constructors

        //- Construct from fvPatch + dgGeomMesh + dictionary of BC file
        dgGeneralBoundaryManager
        (
            const dictionary& fieldDict,
            const dgGeomMesh& mesh
        );

        //- Construct from IOobject + dgGeomMesh
        dgGeneralBoundaryManager
        (
            const IOobject& io,
            const dgGeomMesh& mesh
        );

    //- Destructor
    ~dgGeneralBoundaryManager() = default;

    // Access Functions

        //- Number of boundary patches
        inline label size() const { return bConditions_.size(); }

        //- Const access by patch index
        inline const dgGeneralBoundaryField<Type>& operator[](const label i) const
        {
            return bConditions_[i];
        }

        //- Non-const access
        inline dgGeneralBoundaryField<Type>& operator[](const label i)
        {
            return bConditions_[i];
        }

        //- Return internal field value
        inline const Type& internalValue() const { return internalValue_; }

    // Core methods
        
        //- Update values on plus side according to boundary conditions
        void updateValue
        (
            GaussField<Type>& GField
        );

        //- Update gradients on plus side according to boundary conditions
        void updateGrad
        (
            GaussField<Type>& GField
        );
};

#include "dgGeneralBoundaryManagerI.H"

} // End namespace Foam

#endif
// ************************************************************************* //
