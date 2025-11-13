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

#include "dgGeneralBoundaryManager.H"
#include "tmp.H"
#include "tensor.H"

namespace Foam
{

// * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * * //

template<class Type>
dgGeneralBoundaryManager<Type>::dgGeneralBoundaryManager
(
    const dictionary& fieldDict
)
{
    // Read internalField
    const entry* internalFieldEntry = fieldDict.findEntry("internalField", keyType::LITERAL);

    if (!internalFieldEntry)
    {
        FatalIOErrorInFunction(fieldDict)
            << "Missing 'internalField' entry in field dictionary"
            << exit(FatalIOError);
    }

    // Remove const to allow reading
    ITstream& is = const_cast<entry*>(internalFieldEntry)->stream();

    word fieldType;
    is >> fieldType;

    if (fieldType != "uniform")
    {
        FatalIOErrorInFunction(fieldDict)
            << "Only 'uniform' internalField is supported. Got: " << fieldType
            << exit(FatalIOError);
    }

    is >> internalValue_;

    // Read boundaryField
    const dictionary& bfDict = fieldDict.subDict("boundaryField");
    bConditions_.setSize(bfDict.size());

    label patchI = 0;
    for (const entry& e : bfDict)
    {
        if (!e.isDict()) continue;

        const word& patchName = e.keyword();
        const dictionary& patchDict = e.dict();

        autoPtr<dgGeneralBoundaryField<Type>> bc =
            dgGeneralBoundaryField<Type>::New(patchName, patchDict);

        bConditions_.set(patchI++, bc.ptr());
    }
}

template<class Type>
dgGeneralBoundaryManager<Type>::dgGeneralBoundaryManager
(
    const IOobject& io
)
{
    // Create dictionary from IOobject (e.g., 0/U, 0/p)
    IOdictionary fieldDict(io);

    // Read internalField
    const entry* internalFieldEntry = fieldDict.findEntry("internalField", keyType::LITERAL);

    if (!internalFieldEntry)
    {
        FatalIOErrorInFunction(fieldDict)
            << "Missing 'internalField' entry in field dictionary"
            << exit(FatalIOError);
    }

    // Remove const to allow reading
    ITstream& is = const_cast<entry*>(internalFieldEntry)->stream();

    word fieldType;
    is >> fieldType;

    if (fieldType != "uniform")
    {
        FatalIOErrorInFunction(fieldDict)
            << "Only 'uniform' internalField is supported. Got: " << fieldType
            << exit(FatalIOError);
    }

    is >> internalValue_;

    // Read boundaryField
    const dictionary& bfDict = fieldDict.subDict("boundaryField");
    bConditions_.setSize(bfDict.size());

    label patchI = 0;
    for (const entry& e : bfDict)
    {
        if (!e.isDict()) continue;

        const word& patchName = e.keyword();
        const dictionary& patchDict = e.dict();

        autoPtr<dgGeneralBoundaryField<Type>> bc =
            dgGeneralBoundaryField<Type>::New(patchName, patchDict);

        bConditions_.set(patchI++, bc.ptr());
    }
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

/*
template<class Type>
void dgGeneralBoundaryManager<Type>::setContext(const fieldsContext& ctx)
{
    forAll(bConditions_, i)
    {
        bConditions_[i].setContext(ctx);
    }
}
*/
// * * * * * * * * * * * * * Template Instantiations * * * * * * * * * * * * //

template class dgGeneralBoundaryManager<scalar>;
template class dgGeneralBoundaryManager<vector>;
template class dgGeneralBoundaryManager<tensor>;

// ************************************************************************* //

} // End namespace Foam

