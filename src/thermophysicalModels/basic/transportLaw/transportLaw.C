/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2025
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

#include "transportLaw.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    // Register type name and default debug switch
    defineTypeNameAndDebug(transportLaw, 0);

    // Define runtime selection table for dictionary-based construction
    defineRunTimeSelectionTable(transportLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
 * Constructor:
 * - Store model name and coefficient dictionary snapshot (non-empty by design).
 * - No parameter parsing here; derived classes override read().
 */
Foam::transportLaw::transportLaw
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
:
    name_(name),
    coeff_(dict),
    mesh_(mesh),
    thermo_(thermo)
{
    // Nothing else; derived classes call read() to parse coeff_ fields.
}


// * * * * * * * * * * * * * * * *  Factory Method  * * * * * * * * * * * * * //

/*
 * New(name, dict):
 * - Lookup constructor in the runtime selection table.
 * - If not found, list valid types and abort.
 */
Foam::autoPtr<Foam::transportLaw> Foam::transportLaw::New
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
{
    // Find the matching constructor in the runtime table
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown transportLaw type: " << name << nl
            << "Valid transportLaw types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // Construct and return the selected model
    return cstrIter()(name, dict, mesh, thermo);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

/*
 * read():
 * - Base no-op. Derived classes should:
 *     * extract required parameters from coeff()
 *     * validate ranges/units where appropriate
 */
void Foam::transportLaw::read()
{
    // No-op in base.
}

// ************************************************************************* //

} // End namespace Foam
