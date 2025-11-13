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

#include "dgGeneralBoundaryField.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(dgGeneralBoundaryField<scalar>, 0);
defineNamedTemplateTypeNameAndDebug(dgGeneralBoundaryField<vector>, 0);
defineNamedTemplateTypeNameAndDebug(dgGeneralBoundaryField<tensor>, 0);

defineTemplateRunTimeSelectionTable(dgGeneralBoundaryField<scalar>, dictionary);
defineTemplateRunTimeSelectionTable(dgGeneralBoundaryField<vector>, dictionary);
defineTemplateRunTimeSelectionTable(dgGeneralBoundaryField<tensor>, dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::dgGeneralBoundaryField<Type>::dgGeneralBoundaryField
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    dict_(dict)
{}


// * * * * * * * * * * * * * Factory Method * * * * * * * * * * * * * * //

template<class Type>
autoPtr<dgGeneralBoundaryField<Type>> Foam::dgGeneralBoundaryField<Type>::New
(
    const word& name,
    const dictionary& dict
)
{
    const word bcType(dict.get<word>("type"));
    auto cstrIter = dictionaryConstructorTablePtr_->find(bcType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgGeneralBoundaryField type: " << bcType << nl
            << "Valid types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dict);
}


// * * * * * * * * * * * * * Template Instantiations * * * * * * * * * * * * //

template class dgGeneralBoundaryField<scalar>;
template class dgGeneralBoundaryField<vector>;
template class dgGeneralBoundaryField<tensor>;

} // End namespace Foam

// ************************************************************************* //

