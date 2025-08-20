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
    Foam::dgFixedValueBoundaryField

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
    dgFixedValueBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgFixedValueBoundaryField_H
#define dgFixedValueBoundaryField_H

#include "dgBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgFixedValueBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgFixedValueBoundaryField
:
    public dgBoundaryField<Type>
{
    //- Prescribed boundary value
    Type value_;

public:

    //- Runtime type information
    TypeName("fixedValue");

    //- Constructor from name and dictionary
    dgFixedValueBoundaryField
    (
        const word& name,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgFixedValueBoundaryField() = default;

    //- Update value at plus side using fixed value BC
    virtual void updateValue
    (
        const vector& n,
        const Type& minusValue,
        const Type& minusGrad,
        Type& plusValue,
        Type& plusGrad
    ) const override;

    //- Update gradient at plus side using mirror condition
    virtual void updateGrad
    (
        const vector& n,
        const Type& minusValue,
        const Type& minusGrad,
        Type& plusValue,
        Type& plusGrad
    ) const override;
};

} // End namespace Foam

#endif

