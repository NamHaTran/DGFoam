/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "eqnOfState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // Register type name and default debug switch
    defineTypeNameAndDebug(eqnOfState, 0);

    // Define runtime selection table for dictionary-based construction
    defineRunTimeSelectionTable(eqnOfState, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/* Construct from model name and coefficients dictionary.
 * - Store identifiers; derived classes parse/validate their own keys.
 */
Foam::eqnOfState::eqnOfState
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),   // store model id for logging/debug
    dict_(dict)    // retain coefficients dictionary (const view)
{}


// * * * * * * * * * * * * * * * *  Factory  * * * * * * * * * * * * * * * * //

/* Create EOS by looking up the constructor in the dictionary table.
 * - On failure, list valid types and abort with IO context.
 */
Foam::autoPtr<Foam::eqnOfState> Foam::eqnOfState::New
(
    const word& name,
    const dictionary& dict
)
{
    // 1) Find constructor functor by key (model name)
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    // 2) Report helpful error if not found
    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown eqnOfState type: " << name << nl
            << "Valid eqnOfState types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // 3) Invoke constructor and return autoPtr
    return cstrIter()(name, dict);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

// ************************************************************************* //
