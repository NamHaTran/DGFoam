/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
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

#include "dgCompressibleTemperatureJumpModel.H"
#include "dgThermoConservative.H"
#include "fvPatch.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleTemperatureJumpModel, 0);
defineRunTimeSelectionTable(dgCompressibleTemperatureJumpModel, dictionary);

dgCompressibleTemperatureJumpModel::dgCompressibleTemperatureJumpModel
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    patch_(patch),
    dgMesh_(dgMesh),
    thermo_(thermo),
    dict_(dict)
{}


autoPtr<dgCompressibleTemperatureJumpModel>
dgCompressibleTemperatureJumpModel::New
(
    const word& modelName,
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(modelName);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown temperature jump model: " << modelName << nl
            << "Valid temperature jump models are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(patch, dgMesh, thermo, dict);
}

} // End namespace Foam

// ************************************************************************* //
