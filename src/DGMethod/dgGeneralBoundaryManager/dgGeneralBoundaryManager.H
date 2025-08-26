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

#include "dictionary.H"
#include "wordList.H"
#include "PtrList.H"
#include "dgGeneralBoundaryField.H"
#include "fieldsContext.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgGeneralBoundaryManager Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgGeneralBoundaryManager
{
private:

    //- List of boundary condition objects (one per patch)
    PtrList<dgGeneralBoundaryField<Type>> bConditions_;

    //- Internal field value (only supports 'uniform <value>')
    Type internalValue_;

public:

    // Constructors

        //- Construct from dictionary, without context
        dgGeneralBoundaryManager(const dictionary& fieldDict);

        //- Construct from IOobject
        dgGeneralBoundaryManager(const IOobject& io);

    //- Destructor
    ~dgGeneralBoundaryManager() = default;

    // Context

        //- Set context for all boundary conditions
        void setContext(const fieldsContext& ctx);

    // Access

        //- Return number of patches
        inline label size() const { return bConditions_.size(); }

        //- Access to boundary condition for patch i
        inline const dgGeneralBoundaryField<Type>& operator[](const label i) const
        {
            return bConditions_[i];
        }

        //- Non-const access if needed
        inline dgGeneralBoundaryField<Type>& operator[](const label i)
        {
            return bConditions_[i];
        }

        //- Return internal value
        inline const Type& internalValue() const { return internalValue_; }
};

} // End namespace Foam

//#ifdef NoRepository
//#include "dgGeneralBoundaryManager.C"
//#endif

#endif
