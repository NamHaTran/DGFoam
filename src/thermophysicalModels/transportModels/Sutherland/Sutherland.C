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

#include "Sutherland.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H"

#include <cmath>

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Sutherland, 0);
addToRunTimeSelectionTable(transportLaw, Sutherland, dictionary);


// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

Sutherland::Sutherland
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
:
    transportLaw(name, dict, mesh, thermo),
    As_(1.458e-5),
    S_(110.4),
    Pr0_(0.72)
{
    read();
}


// * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

Foam::scalar Sutherland::calcMu(const scalar T) const
{
    // Guard: T must be positive
    if (T <= scalar(0))
    {
        FatalErrorInFunction
            << "Non-positive temperature in Sutherland::calcMu()"
            << nl << exit(FatalError);
    }

    return As_ * std::pow(T, 1.5)/(T + S_);
}


Foam::scalar Sutherland::calcKappa(const scalar T) const
{
    return calcMu(T)*thermo_.calcCp(T)/Pr0_;
}


Foam::scalar Sutherland::calcPr(const scalar) const
{
    return Pr0_;
}


void Sutherland::read()
{
    if (coeff_.found("As"))
    {
        As_ = readScalar(coeff_.lookup("As"));
    }

    if (coeff_.found("S"))
    {
        S_ = readScalar(coeff_.lookup("S"));
    }

    if (coeff_.found("Pr"))
    {
        Pr0_ = readScalar(coeff_.lookup("Pr"));
    }

    if (As_ <= 0)
    {
        FatalErrorInFunction
            << "Invalid Sutherland 'As'. Must be > 0."
            << exit(FatalError);
    }

    if (S_ <= 0)
    {
        FatalErrorInFunction
            << "Invalid Sutherland 'S'. Must be > 0."
            << exit(FatalError);
    }

    if (Pr0_ <= 0)
    {
        FatalErrorInFunction
            << "Invalid Prandtl number. Must be > 0."
            << exit(FatalError);
    }
}

} // End namespace Foam
