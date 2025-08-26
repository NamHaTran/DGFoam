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
    Foam::dgGeneralZeroGradientBoundaryField

Description
    Implements a zero-gradient (Neumann-type) general boundary condition for 
    Discontinuous Galerkin methods.

    This condition assumes that the field has zero normal derivative across the 
    boundary face. It is implemented by mirroring the internal value to the 
    outside:

        U⁺ = U⁻
        ∇U⁺ = ∇U⁻

    where:
        - U⁻: solution inside the domain
        - U⁺: ghost value used for flux computation
        - ∇U⁻: gradient inside the domain
        - ∇U⁺: gradient used for flux evaluation on ghost side

    Typical use:
        - Outlet boundaries
        - Symmetry planes (for scalar fields)
        - Free-slip walls (when velocity gradient is zero)

SourceFiles
    dgGeneralZeroGradientBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralZeroGradientBoundaryField_H
#define dgGeneralZeroGradientBoundaryField_H

#include "dgGeneralBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class dgGeneralZeroGradientBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgGeneralZeroGradientBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
public:

    //- Runtime type information
    TypeName("generalZeroGradient");

    //- Constructor from name and dictionary
    dgGeneralZeroGradientBoundaryField
    (
        const word& name,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralZeroGradientBoundaryField() = default;

    //- Update plus-side value using zero-gradient BC: U+ = U-, S+ = S-
    virtual void updateValue
    (
        const label gaussID,
        const vector& n,
        const Type& minusValue,
        const Type& minusGrad,
        Type& plusValue,
        Type& plusGrad
    ) const override;

    //- Update plus-side gradient (copy from minus side)
    virtual void updateGrad
    (
        const label gaussID,
        const vector& n,
        const Type& minusValue,
        const Type& minusGrad,
        Type& plusValue,
        Type& plusGrad
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //

