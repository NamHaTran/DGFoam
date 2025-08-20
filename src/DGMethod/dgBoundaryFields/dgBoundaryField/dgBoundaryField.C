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

#include "dgBoundaryField.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Type info for common types
defineNamedTemplateTypeNameAndDebug(dgBoundaryField<scalar>, 0);
defineNamedTemplateTypeNameAndDebug(dgBoundaryField<vector>, 0);
defineNamedTemplateTypeNameAndDebug(dgBoundaryField<tensor>, 0);

// Template specializations for runtime selection table
defineTemplateRunTimeSelectionTable(dgBoundaryField<scalar>, dictionary);
defineTemplateRunTimeSelectionTable(dgBoundaryField<vector>, dictionary);
defineTemplateRunTimeSelectionTable(dgBoundaryField<tensor>, dictionary);


// * * * * * * * * * * * * * Private/Helper Functions  * * * * * * * * * * * //

template<class Type>
dgBCForm Foam::dgBoundaryField<Type>::parseForm_(const dictionary& dict)
{
    // 1) Read 'form' from dict with default "strong"
    const word f = dict.lookupOrDefault<word>("form", "strong");

    // 2) Map to enum; accept case-insensitive tokens
    if (f == "strong" || f == "Strong" || f == "STRONG")
    {
        return dgBCForm::strong;
    }
    else if (f == "weak" || f == "Weak" || f == "WEAK")
    {
        return dgBCForm::weak;
    }

    // 3) Fallback + diagnostic
    WarningInFunction
        << "Unknown form '" << f << "'. Using 'strong' by default.\n";

    return dgBCForm::strong;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::dgBoundaryField<Type>::dgBoundaryField
(
    const word& name,
    const dictionary& dict
)
:
    // 1) Copy basic inputs
    name_(name),
    dict_(dict),
    ctxPtr_(nullptr),

    // 2) Parse and store enforcement form
    bcForm_(parseForm_(dict_))
{
    // Note:
    // - ctxPtr_ is set later via setContext() when available.
    // - bcForm_ determines how derived classes implement updateValue/updateGrad.
}


// * * * * * * * * * * * * * * * *  Factory Method  * * * * * * * * * * * * * //

template<class Type>
autoPtr<dgBoundaryField<Type>> Foam::dgBoundaryField<Type>::New
(
    const word& name,
    const dictionary& dict
)
{
    // 1) Get 'type' and find constructor
    const word bcType(dict.get<word>("type"));
    auto cstrIter = dictionaryConstructorTablePtr_->find(bcType);

    // 2) Validate
    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgBoundaryField type: " << bcType << nl
            << "Valid types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // 3) Construct and return the selected model
    return cstrIter()(name, dict);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class Type>
word Foam::dgBoundaryField<Type>::formName() const
{
    // Return a stable token for the current form
    return (bcForm_ == dgBCForm::strong ? word("strong") : word("weak"));
}


// * * * * * * * * * * * * * * Template instantiation * * * * * * * * * * * * //

template class dgBoundaryField<scalar>;
template class dgBoundaryField<vector>;
template class dgBoundaryField<tensor>;

template autoPtr<dgBoundaryField<scalar>>
dgBoundaryField<scalar>::New(const word&, const dictionary&);

template autoPtr<dgBoundaryField<vector>>
dgBoundaryField<vector>::New(const word&, const dictionary&);

template autoPtr<dgBoundaryField<tensor>>
dgBoundaryField<tensor>::New(const word&, const dictionary&);

} // End namespace Foam

// ************************************************************************* //
