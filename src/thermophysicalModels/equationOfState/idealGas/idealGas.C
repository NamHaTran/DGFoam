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

#include "idealGas.H"
#include "addToRunTimeSelectionTable.H" 
#include "constants.H"

namespace Foam
{
    // Register type name and add to runtime selection table
    defineTypeNameAndDebug(idealGas, 0);
    addToRunTimeSelectionTable(eqnOfState, idealGas, dictionary);
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

Foam::idealGas::idealGas
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    eqnOfState(name, dict, mesh),
    R_
    (
        constant::physicoChemical::R.value()
      / readScalar(dict.subDict("specie").lookup("molWeight"))
    )
{}


// * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

// rho = p / (R * T)
void Foam::idealGas::calcRhoFromPT
(
    const label cellI,
    const GaussField<scalar>& p,
    const GaussField<scalar>& T,
    GaussField<scalar>& rho
) const
{
    rho = p / (R_ * T);
}


// p = rho * R * T
void Foam::idealGas::calcPFromRhoT
(
    const label cellI,
    const GaussField<scalar>& rho,
    const GaussField<scalar>& T,
    GaussField<scalar>& p
) const
{
    p = rho * (R_ * T);
}


// T = p / (rho * R)
void Foam::idealGas::calcTFromPRho
(
    const label cellI,
    const GaussField<scalar>& p,
    const GaussField<scalar>& rho,
    GaussField<scalar>& T
) const
{
    T = p / (rho * R_);
}

// ************************************************************************* //
