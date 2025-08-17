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

#include "thermoLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // Register type name and default debug switch
    defineTypeNameAndDebug(thermoLaw, 0);

    // Define runtime selection table for dictionary-based construction
    defineRunTimeSelectionTable(thermoLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
 * Constructor: stores the model name and coefficient dictionary.
 * Any parameter extraction/validation is deferred to the derived class.
 */
Foam::thermoLaw::thermoLaw
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * *  Factory Method  * * * * * * * * * * * * * //

/*
 * New(name, dict):
 * - Look up the constructor from the runtime selection table.
 * - Create and return a new thermoLaw instance.
 * - If not found, print available types and exit with an error.
 */
Foam::autoPtr<Foam::thermoLaw> Foam::thermoLaw::New
(
    const word& name,
    const dictionary& dict
)
{
    // Find the matching constructor in the runtime table
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown thermoLaw type: " << name << nl
            << "Valid thermoLaw types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // Construct and return the selected model
    return cstrIter()(name, dict);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::thermoLaw::a(const scalar) const
{
    FatalErrorInFunction
        << "Speed of sound a(T) not implemented for thermoLaw type: " << name_
        << exit(FatalError);

    // Unreachable; return to satisfy compilers
    return 0.0;
}

// ************************************************************************* //
