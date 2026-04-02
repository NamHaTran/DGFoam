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

#include "thermoLaw.H"
#include "dgExpr.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // Register type name and default debug switch
    defineTypeNameAndDebug(thermoLaw, 0);

    // Define runtime selection table for dictionary-based construction
    defineRunTimeSelectionTable(thermoLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
 * Constructor: stores the model name and coefficient dictionary.
 * Any parameter extraction/validation is deferred to the derived class.
 */
Foam::thermoLaw::thermoLaw
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const eqnOfState& eos
)
:
    name_(name),
    dict_(dict),
    mesh_(mesh),
    eos_(eos)
{}


// * * * * * * * * * * * * * * * *  Factory Method  * * * * * * * * * * * * * //

/*
 * New(name, dict):
 * - Look up the constructor from the runtime selection table.
 * - Create and return a new thermoLaw instance.
 * - If not found, print available types and exit with an error.
 */
Foam::autoPtr<Foam::thermoLaw> Foam::thermoLaw::New
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const eqnOfState& eos
)
{
    // Find the matching constructor in the runtime table
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown thermoLaw type: " << name << nl
            << "Valid thermoLaw types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // Construct and return the selected model
    return cstrIter()(name, dict, mesh, eos);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::thermoLaw::calcTFromInternalE
(
    const label,
    const GaussField<scalar>& e,
    GaussField<scalar>& T
) const
{
    dg::assign
    (
        T,
        dg::map
        (
            [this](const scalar eValue)
            {
                return calcTFromInternalE(eValue);
            },
            dg::expr(e)
        )
    );
}


void Foam::thermoLaw::calcTFromInternalE
(
    const boundaryGaussField<scalar>& e,
    boundaryGaussField<scalar>& T
) const
{
    for (label gpI = 0; gpI < e.size(); ++gpI)
    {
        T[gpI] = calcTFromInternalE(e[gpI]);
    }
}


void Foam::thermoLaw::calcTFromH
(
    const label,
    const GaussField<scalar>& h,
    GaussField<scalar>& T
) const
{
    dg::assign
    (
        T,
        dg::map
        (
            [this](const scalar hValue)
            {
                return calcTFromH(hValue);
            },
            dg::expr(h)
        )
    );
}


void Foam::thermoLaw::calcTFromH
(
    const boundaryGaussField<scalar>& h,
    boundaryGaussField<scalar>& T
) const
{
    for (label gpI = 0; gpI < h.size(); ++gpI)
    {
        T[gpI] = calcTFromH(h[gpI]);
    }
}


void Foam::thermoLaw::calcCp
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& Cp
) const
{
    dg::assign
    (
        Cp,
        dg::map
        (
            [this](const scalar TValue)
            {
                return calcCp(TValue);
            },
            dg::expr(T)
        )
    );
}


void Foam::thermoLaw::calcCp
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& Cp
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        Cp[gpI] = calcCp(T[gpI]);
    }
}


void Foam::thermoLaw::calcCv
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& Cv
) const
{
    dg::assign
    (
        Cv,
        dg::map
        (
            [this](const scalar TValue)
            {
                return calcCv(TValue);
            },
            dg::expr(T)
        )
    );
}


void Foam::thermoLaw::calcCv
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& Cv
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        Cv[gpI] = calcCv(T[gpI]);
    }
}


void Foam::thermoLaw::calcGamma
(
    const label,
    const GaussField<scalar>& Cp,
    const GaussField<scalar>& Cv,
    GaussField<scalar>& gamma
) const
{
    dg::assign
    (
        gamma,
        dg::map
        (
            [this](const scalar cpValue, const scalar cvValue)
            {
                return calcGamma(cpValue, cvValue);
            },
            dg::expr(Cp),
            dg::expr(Cv)
        )
    );
}


void Foam::thermoLaw::calcGamma
(
    const boundaryGaussField<scalar>& Cp,
    const boundaryGaussField<scalar>& Cv,
    boundaryGaussField<scalar>& gamma
) const
{
    for (label gpI = 0; gpI < Cp.size(); ++gpI)
    {
        gamma[gpI] = calcGamma(Cp[gpI], Cv[gpI]);
    }
}


void Foam::thermoLaw::calcInternalE
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& e
) const
{
    dg::assign
    (
        e,
        dg::map
        (
            [this](const scalar TValue)
            {
                return calcInternalE(TValue);
            },
            dg::expr(T)
        )
    );
}


void Foam::thermoLaw::calcInternalE
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& e
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        e[gpI] = calcInternalE(T[gpI]);
    }
}


void Foam::thermoLaw::calcH
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& h
) const
{
    dg::assign
    (
        h,
        dg::map
        (
            [this](const scalar TValue)
            {
                return calcH(TValue);
            },
            dg::expr(T)
        )
    );
}


void Foam::thermoLaw::calcH
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& h
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        h[gpI] = calcH(T[gpI]);
    }
}


void Foam::thermoLaw::calcSpeedOfSound
(
    const label,
    const GaussField<scalar>& T,
    const GaussField<scalar>& gamma,
    GaussField<scalar>& a
) const
{
    dg::assign
    (
        a,
        dg::map
        (
            [this](const scalar TValue, const scalar gammaValue)
            {
                return calcSpeedOfSound(TValue, gammaValue);
            },
            dg::expr(T),
            dg::expr(gamma)
        )
    );
}


void Foam::thermoLaw::calcSpeedOfSound
(
    const boundaryGaussField<scalar>& T,
    const boundaryGaussField<scalar>& gamma,
    boundaryGaussField<scalar>& a
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        a[gpI] = calcSpeedOfSound(T[gpI], gamma[gpI]);
    }
}


// ************************************************************************* //
