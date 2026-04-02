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

#include "eqnOfState.H"
#include "dgExpr.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // Register type name and default debug switch
    defineTypeNameAndDebug(eqnOfState, 0);

    // Define runtime selection table for dictionary-based construction
    defineRunTimeSelectionTable(eqnOfState, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/* Construct from model name and coefficients dictionary.
 * - Store identifiers; derived classes parse/validate their own keys.
 */
Foam::eqnOfState::eqnOfState
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    name_(name),   // store model id for logging/debug
    dict_(dict),   // retain coefficients dictionary (const view)
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * *  Factory  * * * * * * * * * * * * * * * * //

/* Create EOS by looking up the constructor in the dictionary table.
 * - On failure, list valid types and abort with IO context.
 */
Foam::autoPtr<Foam::eqnOfState> Foam::eqnOfState::New
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    // 1) Find constructor functor by key (model name)
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    // 2) Report helpful error if not found
    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown eqnOfState type: " << name << nl
            << "Valid eqnOfState types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // 3) Invoke constructor and return autoPtr
    return cstrIter()(name, dict, mesh);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::eqnOfState::calcRhoFromPT
(
    const label,
    const GaussField<scalar>& p,
    const GaussField<scalar>& T,
    GaussField<scalar>& rho
) const
{
    dg::assign
    (
        rho,
        dg::map
        (
            [this](const scalar pValue, const scalar TValue)
            {
                return calcRhoFromPT(pValue, TValue);
            },
            dg::expr(p),
            dg::expr(T)
        )
    );
}


void Foam::eqnOfState::calcRhoFromPT
(
    const boundaryGaussField<scalar>& p,
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& rho
) const
{
    for (label gpI = 0; gpI < p.size(); ++gpI)
    {
        rho[gpI] = calcRhoFromPT(p[gpI], T[gpI]);
    }
}


void Foam::eqnOfState::calcPFromRhoT
(
    const label,
    const GaussField<scalar>& rho,
    const GaussField<scalar>& T,
    GaussField<scalar>& p
) const
{
    dg::assign
    (
        p,
        dg::map
        (
            [this](const scalar rhoValue, const scalar TValue)
            {
                return calcPFromRhoT(rhoValue, TValue);
            },
            dg::expr(rho),
            dg::expr(T)
        )
    );
}


void Foam::eqnOfState::calcPFromRhoT
(
    const boundaryGaussField<scalar>& rho,
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& p
) const
{
    for (label gpI = 0; gpI < rho.size(); ++gpI)
    {
        p[gpI] = calcPFromRhoT(rho[gpI], T[gpI]);
    }
}


void Foam::eqnOfState::calcTFromPRho
(
    const label,
    const GaussField<scalar>& p,
    const GaussField<scalar>& rho,
    GaussField<scalar>& T
) const
{
    dg::assign
    (
        T,
        dg::map
        (
            [this](const scalar pValue, const scalar rhoValue)
            {
                return calcTFromPRho(pValue, rhoValue);
            },
            dg::expr(p),
            dg::expr(rho)
        )
    );
}


void Foam::eqnOfState::calcTFromPRho
(
    const boundaryGaussField<scalar>& p,
    const boundaryGaussField<scalar>& rho,
    boundaryGaussField<scalar>& T
) const
{
    for (label gpI = 0; gpI < p.size(); ++gpI)
    {
        T[gpI] = calcTFromPRho(p[gpI], rho[gpI]);
    }
}

// ************************************************************************* //
