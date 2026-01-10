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
    Foam::dgProcessorBoundaryField

Description


SourceFiles
    dgProcessorBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgProcessorBoundaryField_H
#define dgProcessorBoundaryField_H

#include "dgInterfaceBoundaryField.H"
#include "processorFvPatch.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgProcessorBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgProcessorBoundaryField
:
    public dgInterfaceBoundaryField<Type>
{
private:

    //- Reference to processor patch (type-safe)
    const processorFvPatch& procPatch_;

    //- Neighbour processor rank
    const label neighbProc_;

    //- Communication buffers
    mutable List<Type> sendBuf_;
    mutable List<Type> recvBuf_;

public:

    //- Runtime type information
    TypeName("processor");

    //- Constructor from name and dictionary
    dgProcessorBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh
    );

    //- Destructor
    virtual ~dgProcessorBoundaryField() = default;

    //- Exchange Gauss data across processor patch
    void update
    (
        dgField<Type>& field
    ) const override;

    //- Check patch type
    void checkPatchType() const override
    {
        if (
            this->patch_.type() != "processor"
        )
        {
            FatalErrorInFunction
                << "Boundary condition " << type() <<" can only "
                << "be applied to patch type: \n"
                << "    " << "processor \n"
                << "but patch "
                << this->patch_.name() << " is of type " << this->patch_.type() << nl
                << exit(FatalError);
        }
    }
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgProcessorBoundaryFieldI.H"

#endif
// ************************************************************************* //
