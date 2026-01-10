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
    Foam::dgProcessorBoundaryManager

Description
    Class for DG processor boundary conditions management.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgProcessorBoundaryManager_H
#define Foam_dgProcessorBoundaryManager_H

#include "dgGeomMesh.H"
#include "fvPatch.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "dgProcessorBoundaryField.H"
#include "dgField.H"
#include "dgBoundaryHelper.H"
#include "Pstream.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgProcessorBoundaryManager Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgProcessorBoundaryManager
{
private:

    //- Reference to dgGeomMesh (used for fvPatch access and field lookup)
    const dgGeomMesh& mesh_;

    //- List of boundary condition objects (one per fvPatch)
    PtrList<dgProcessorBoundaryField<Type>> bConditions_;

public:

    // Constructors

        //- Construct from fvPatch + dgGeomMesh
        dgProcessorBoundaryManager
        (
            const dgGeomMesh& mesh
        );

    //- Destructor
    ~dgProcessorBoundaryManager() = default;

    // Access Functions

        //- Number of boundary patches
        inline label size() const { return bConditions_.size(); }

        //- Const access by patch index
        inline const dgProcessorBoundaryField<Type>& operator[](const label i) const
        {
            return bConditions_[i];
        }

        //- Non-const access
        inline dgProcessorBoundaryField<Type>& operator[](const label i)
        {
            return bConditions_[i];
        }

    // Core methods

        //- Update field values at the boundary
        void update
        (
            dgField<Type>& field
        );
};

#include "dgProcessorBoundaryManagerI.H"

} // End namespace Foam

#endif
// ************************************************************************* //
