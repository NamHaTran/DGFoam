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

#include "constantCp.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "IOobject.H"
#include "mathematicalConstants.H"
#include "constants.H"   // for constant::physicoChemical::R

#include <cmath>

namespace Foam
{
    // Register into thermoLaw (dictionary-based construction)
    defineTypeNameAndDebug(constantCp, 0);
    addToRunTimeSelectionTable(thermoLaw, constantCp, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
 * Construct from root dictionary:
 * - Reads Cp from 'thermodynamics' sub-dictionary
 * - Reads molWeight from 'specie' sub-dictionary
 * - Computes R = Ru / molWeight (OF convention: molWeight in kg/kmol,
 *   Ru in J/(kmol.K)), Cv = Cp - R, and gamma = Cp/Cv
 */
Foam::constantCp::constantCp
(
    const word& name,
    const dictionary& rootDict
)
:
    thermoLaw(name, rootDict),
    Cp_(Zero),
    R_(Zero),
    Cv_(Zero),
    gamma_(Zero)
{
    // 1) Read Cp from 'thermodynamics' sub-dictionary
    if (!rootDict.found("thermodynamics"))
    {
        FatalIOErrorInFunction(rootDict)
            << "Missing 'thermodynamics' sub-dictionary for thermo model '"
            << name << "'. Expecting 'thermodynamics { Cp <value>; }'."
            << exit(FatalIOError);
    }

    const dictionary& thermoDict = rootDict.subDict("thermodynamics");

    if (!thermoDict.found("Cp"))
    {
        FatalIOErrorInFunction(thermoDict)
            << "Entry 'Cp' not found in 'thermodynamics' for thermo model '"
            << name << "'."
            << exit(FatalIOError);
    }

    Cp_ = readScalar(thermoDict.lookup("Cp"));

    // 2) Read molWeight from 'specie' sub-dictionary to compute R
    if (!rootDict.found("specie"))
    {
        FatalIOErrorInFunction(rootDict)
            << "Missing 'specie' sub-dictionary. Required to compute R from molWeight."
            << exit(FatalIOError);
    }

    const dictionary& specieDict = rootDict.subDict("specie");

    if (!specieDict.found("molWeight"))
    {
        FatalIOErrorInFunction(specieDict)
            << "Entry 'molWeight' not found in 'specie'. Required to compute R."
            << exit(FatalIOError);
    }

    const scalar W = readScalar(specieDict.lookup("molWeight")); // [kg/kmol] in OF convention

    // 3) Compute R = Ru / W
    //    Note: OpenFOAM constant::physicoChemical::R is in [J/(kmol.K)]
    const scalar Ru = constant::physicoChemical::R.value(); // ~8.314462618e3 J/(kmol.K)
    if (W <= SMALL)
    {
        FatalIOErrorInFunction(specieDict)
            << "Invalid 'molWeight' <= 0."
            << exit(FatalIOError);
    }

    R_ = Ru / W; // [J/(kg.K)]

    // 4) Compute Cv and gamma
    Cv_ = Cp_ - R_;
    if (Cv_ <= SMALL)
    {
        FatalIOErrorInFunction(thermoDict)
            << "Computed Cv = Cp - R <= 0 (Cp=" << Cp_ << ", R=" << R_
            << "). Check your inputs."
            << exit(FatalIOError);
    }

    gamma_ = Cp_ / Cv_;
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::constantCp::Cp(const scalar T) const
{
    // Constant Cp model
    return Cp_;
}

Foam::scalar Foam::constantCp::h(const scalar T) const
{
    // h(T) = Cp * T  (no reference offset here)
    return Cp_ * T;
}

Foam::scalar Foam::constantCp::e(const scalar T) const
{
    // e(T) = Cv * T = (Cp - R) * T
    return Cv_ * T;
}

Foam::scalar Foam::constantCp::TfrH(const scalar h) const
{
    return h / Cp_;
}

Foam::scalar Foam::constantCp::TfrE(const scalar e) const
{
    return e / Cv_;
}

Foam::scalar Foam::constantCp::a(const scalar T) const
{
    return std::sqrt(max(gamma_*R_*T, SMALL));
}
// ************************************************************************* //
