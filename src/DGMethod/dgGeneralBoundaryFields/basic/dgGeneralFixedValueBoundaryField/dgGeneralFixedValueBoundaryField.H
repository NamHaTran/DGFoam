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
    Foam::dgGeneralFixedValueBoundaryField

Description
    Fixed value (Dirichlet) boundary condition for DG methods.

    Applies a strong Dirichlet constraint at the boundary by mirroring
    the internal value around the prescribed boundary value:

        U⁺ = 2·U_bc − U⁻

    The gradient is not modified (copied as-is), assuming the Dirichlet
    value controls the field, not the gradient.

    Supports any Type (scalar, vector, tensor, etc.) that provides basic
    arithmetic operations.

SourceFiles
    dgGeneralFixedValueBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralFixedValueBoundaryField_H
#define dgGeneralFixedValueBoundaryField_H

#include "dgGeneralBoundaryField.H"
#include "wallPolyPatch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgGeneralFixedValueBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgGeneralFixedValueBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
    //- Prescribed boundary value
    Type value_;

public:

    //- Runtime type information
    TypeName("generalFixedValue");

    //- Constructor from name and dictionary
    dgGeneralFixedValueBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralFixedValueBoundaryField() = default;

    //- Update value at plus side using fixed value BC
    virtual void updateGhostStateValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const Type& minusValue,
        Type& plusValue
    ) const override
    {
        plusValue = 2 * value_ - minusValue;
    }

    //- Check patch type
    void checkPatchType() const override
    {
        if (
            this->patch_.type() != "patch" &&
            this->patch_.type() != "wall"
        )
        {
            FatalErrorInFunction
                << "Boundary condition " << type() <<" can only "
                << "be applied to patch type: \n"
                << "    " << "patch \n"
                << "    " << "wall \n"
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

#include "dgGeneralFixedValueBoundaryFieldI.H"

#endif
// ************************************************************************* //