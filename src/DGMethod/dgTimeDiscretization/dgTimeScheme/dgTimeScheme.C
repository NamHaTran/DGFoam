/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgTimeScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{

defineTypeNameAndDebug(dgTimeScheme, 0);
defineRunTimeSelectionTable(dgTimeScheme, dictionary);

// Base wrapper only stores shared metadata.

Foam::dgTimeScheme::dgTimeScheme
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    name_(name),
    dict_(dict),
    mesh_(mesh),
    stageI_(0)
{}


autoPtr<dgTimeScheme> Foam::dgTimeScheme::New
(
    const word& name,
    const word& timeSchemeType,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    // Select the concrete scheme from the runtime table.
    auto cstrIter = dictionaryConstructorTablePtr_->find(timeSchemeType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgTimeScheme type: " << timeSchemeType << nl
            << "Valid types are: " << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dict, mesh);
}

} // End namespace Foam

// ************************************************************************* //
