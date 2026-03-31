/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
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

#include "dgCompressibleBoundaryField.H"
#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"
#include <string>

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleBoundaryField, 0);
defineRunTimeSelectionTable(dgCompressibleBoundaryField, dictionary);

namespace
{

word normalizedCompressibleBoundaryType(const dictionary& dict)
{
    const string rawType(dict.get<string>("type"));
    static const std::string prefix("compressible::");

    if (rawType.size() > prefix.size() && rawType.substr(0, prefix.size()) == prefix)
    {
        return word(rawType);
    }

    return word(prefix + rawType);
}

} // End anonymous namespace

dgCompressibleBoundaryField::dgCompressibleBoundaryField
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


autoPtr<dgCompressibleBoundaryField> dgCompressibleBoundaryField::New
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
{
    const word bcType(normalizedCompressibleBoundaryType(dict));

    auto cstrIter = dictionaryConstructorTablePtr_->find(bcType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgCompressibleBoundaryField type: " << bcType << nl
            << "Valid dgCompressibleBoundaryField types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(patch, dgMesh, thermo, dict);
}


void dgCompressibleBoundaryField::primitiveToConservative
(
    const scalar p,
    const scalar T,
    const vector& U,
    scalar& rho,
    vector& rhoU,
    scalar& E
) const
{
    rho  = thermo_.eos().calcRhoFromPT(p, T);
    rhoU = rho*U;

    const scalar he = thermo_.energyModel().calcHe(T);

    if (thermo_.heIsInternalEnergy())
    {
        E = rho*(he + 0.5*magSqr(U));
    }
    else
    {
        E = rho*(he + 0.5*magSqr(U)) - p;
    }
}


void dgCompressibleBoundaryField::pressureMachTemperatureToConservative
(
    const scalar p,
    const scalar T,
    const vector& Mach,
    scalar& rho,
    vector& rhoU,
    scalar& E
) const
{
    const scalar Cp = thermo_.thermo().calcCp(T);
    const scalar Cv = thermo_.thermo().calcCv(T);
    const scalar gamma = thermo_.thermo().calcGamma(Cp, Cv);
    const scalar a = thermo_.thermo().calcSpeedOfSound(T, gamma);

    primitiveToConservative(p, T, a*Mach, rho, rhoU, E);
}

} // End namespace Foam

// ************************************************************************* //
