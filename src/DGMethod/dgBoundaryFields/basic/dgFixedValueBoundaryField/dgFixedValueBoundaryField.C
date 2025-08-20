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

\*---------------------------------------------------------------------------*/

#include "dgFixedValueBoundaryField.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
// Type info for common types
defineNamedTemplateTypeNameAndDebug(dgFixedValueBoundaryField<scalar>, 0);
defineNamedTemplateTypeNameAndDebug(dgFixedValueBoundaryField<vector>, 0);
defineNamedTemplateTypeNameAndDebug(dgFixedValueBoundaryField<tensor>, 0);

// Register the dgFixedValueBoundaryField into the runtime selection table
addToRunTimeSelectionTable(dgBoundaryField<scalar>, dgFixedValueBoundaryField<scalar>, dictionary);
addToRunTimeSelectionTable(dgBoundaryField<vector>, dgFixedValueBoundaryField<vector>, dictionary);
addToRunTimeSelectionTable(dgBoundaryField<tensor>, dgFixedValueBoundaryField<tensor>, dictionary);


// * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

template<class Type>
dgFixedValueBoundaryField<Type>::dgFixedValueBoundaryField
(
    const word& name,
    const dictionary& dict
)
:
    dgBoundaryField<Type>(name, dict),
    value_(dict.lookup("value"))
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class Type>
void dgFixedValueBoundaryField<Type>::updateValue
(
    const vector& n,
    const Type& minusValue,
    const Type& minusGrad,
    Type& plusValue,
    Type& plusGrad
) const
{
    // Dirichlet condition: mirror around boundary value
    plusValue = 2 * value_ - minusValue;
}

template<class Type>
void dgFixedValueBoundaryField<Type>::updateGrad
(
    const vector& n,
    const Type& minusValue,
    const Type& minusGrad,
    Type& plusValue,
    Type& plusGrad
) const
{
    // Gradient is unchanged
    plusGrad = minusGrad;
}

// * * * * * * * * * * * * * * Template instantiation * * * * * * * * * * * * //
// not necessary

} // End namespace Foam

// ************************************************************************* //
