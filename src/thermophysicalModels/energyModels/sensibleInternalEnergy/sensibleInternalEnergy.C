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

// * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensibleInternalEnergy, 0);
addToRunTimeSelectionTable(energy, sensibleInternalEnergy, dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

sensibleInternalEnergy::sensibleInternalEnergy
(
    const word& name,
    const dictionary& dict
)
:
    energy(name, dict),
    coeffDict_(dict)  // Store in case future coeffs are needed
{}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

scalar sensibleInternalEnergy::T
(
    const scalar he,
    const scalar Cv,
    const scalar Cp
) const
{
    return he / Cv;
}

} // namespace Foam


