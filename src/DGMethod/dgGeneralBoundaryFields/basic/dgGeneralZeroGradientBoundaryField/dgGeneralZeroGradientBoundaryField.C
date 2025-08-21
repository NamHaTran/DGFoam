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

#include "dgGeneralZeroGradientBoundaryField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * Runtime Type * * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(dgGeneralZeroGradientBoundaryField<scalar>, 0);
defineNamedTemplateTypeNameAndDebug(dgGeneralZeroGradientBoundaryField<vector>, 0);
defineNamedTemplateTypeNameAndDebug(dgGeneralZeroGradientBoundaryField<tensor>, 0);

addTemplatedToRunTimeSelectionTable
(
    dgGeneralBoundaryField,
    dgGeneralZeroGradientBoundaryField,
    scalar,
    dictionary
);

addTemplatedToRunTimeSelectionTable
(
    dgGeneralBoundaryField,
    dgGeneralZeroGradientBoundaryField,
    vector,
    dictionary
);

addTemplatedToRunTimeSelectionTable
(
    dgGeneralBoundaryField,
    dgGeneralZeroGradientBoundaryField,
    tensor,
    dictionary
);


// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
dgGeneralZeroGradientBoundaryField<Type>::dgGeneralZeroGradientBoundaryField
(
    const word& name,
    const dictionary& dict
)
:
    dgGeneralBoundaryField<Type>(name, dict)
{}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class Type>
void dgGeneralZeroGradientBoundaryField<Type>::updateValue
(
    const vector& n,
    const Type& minusValue,
    const Type& minusGrad,
    Type& plusValue,
    Type& plusGrad
) const
{
    plusValue = minusValue;
    // gradient is left untouched
}


template<class Type>
void dgGeneralZeroGradientBoundaryField<Type>::updateGrad
(
    const vector& n,
    const Type& minusValue,
    const Type& minusGrad,
    Type& plusValue,
    Type& plusGrad
) const
{
    plusGrad = minusGrad;
    // value is left untouched
}


// * * * * * * * * * * * * * Template Instantiation * * * * * * * * * * * * //

} // End namespace Foam

