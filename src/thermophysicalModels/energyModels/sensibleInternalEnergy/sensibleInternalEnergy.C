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
    energy(name, dict, mesh, thermo),
    coeffDict_(dict)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sensibleInternalEnergy::calcEnthalpy
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& h
) const
{
    // Use thermoLaw to compute:
    //     h = h(T)
    thermo_.calcH(cellI, T, h);
}


void Foam::sensibleInternalEnergy::calcEnergy
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& e
) const
{
    // Use thermoLaw to compute:
    //     e = e(T)
    thermo_.calcInternalE(cellI, T, e);
}


void Foam::sensibleInternalEnergy::calcTfromEnergy
(
    const label cellI,
    const GaussField<scalar>& e,
    GaussField<scalar>& T
) const
{
    // Use thermoLaw to invert:
    //     T = T(e)
    thermo_.calcT(cellI, e, T);
}


// ************************************************************************* //

} // End namespace Foam

