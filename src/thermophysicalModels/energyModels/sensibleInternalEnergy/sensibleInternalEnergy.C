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

#include "sensibleInternalEnergy.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensibleInternalEnergy, 0);
addToRunTimeSelectionTable(energy, sensibleInternalEnergy, dictionary);


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::sensibleInternalEnergy::sensibleInternalEnergy
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
:
    energy(name, dict, mesh, thermo, energy::heType::internalEnergy),
    coeffDict_(dict)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::sensibleInternalEnergy::calcEnthalpy(const scalar T) const
{
    return thermo_.calcH(T);
}


Foam::scalar Foam::sensibleInternalEnergy::calcHe(const scalar T) const
{
    return thermo_.calcInternalE(T);
}


Foam::scalar Foam::sensibleInternalEnergy::calcTfromHe(const scalar he) const
{
    return thermo_.calcTFromInternalE(he);
}


Foam::vector Foam::sensibleInternalEnergy::calcGradTfromHe
(
    const scalar he,
    const vector& gradHe
) const
{
    const scalar T = calcTfromHe(he);
    return gradHe/max(thermo_.calcCv(T), SMALL);
}


// ************************************************************************* //

} // End namespace Foam
