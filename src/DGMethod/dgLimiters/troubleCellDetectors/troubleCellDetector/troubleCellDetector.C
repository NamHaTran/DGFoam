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

#include "troubleCellDetector.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{

defineTypeNameAndDebug(troubleCellDetector, 0);
defineRunTimeSelectionTable(troubleCellDetector, dictionary);


word troubleCellDetector::lookupFieldName
(
    const dictionary& dict,
    const word& entryName,
    const word& defaultName
)
{
    if (dict.found(entryName))
    {
        return dict.get<word>(entryName);
    }

    return defaultName;
}


troubleCellDetector::troubleCellDetector
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    rho_
    (
        mesh_.getFvMesh().lookupObject<dgField<scalar>>
        (
            lookupFieldName(dict, "rho", "rho")
        )
    ),
    rhoU_
    (
        mesh_.getFvMesh().lookupObject<dgField<vector>>
        (
            lookupFieldName(dict, "rhoU", "rhoU")
        )
    ),
    E_
    (
        mesh_.getFvMesh().lookupObject<dgField<scalar>>
        (
            lookupFieldName(dict, "E", "E")
        )
    )
{}


autoPtr<troubleCellDetector> troubleCellDetector::New
(
    const word& detectorType,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(detectorType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown troubleCellDetector type: " << detectorType << nl
            << "Valid troubleCellDetector types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(dict, mesh);
}

} // End namespace Foam

// ************************************************************************* //
