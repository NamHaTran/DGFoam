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

#include "transportLaw.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    // Register type name and default debug switch
    defineTypeNameAndDebug(transportLaw, 0);

    // Define runtime selection table for dictionary-based construction
    defineRunTimeSelectionTable(transportLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
 * Constructor:
 * - Store model name and coefficient dictionary snapshot (non-empty by design).
 * - No parameter parsing here; derived classes override read().
 */
Foam::transportLaw::transportLaw
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
:
    name_(name),
    coeff_(dict),
    mesh_(mesh),
    thermo_(thermo)
{
    // Nothing else; derived classes call read() to parse coeff_ fields.
}


// * * * * * * * * * * * * * * * *  Factory Method  * * * * * * * * * * * * * //

/*
 * New(name, dict):
 * - Lookup constructor in the runtime selection table.
 * - If not found, list valid types and abort.
 */
Foam::autoPtr<Foam::transportLaw> Foam::transportLaw::New
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
{
    // Find the matching constructor in the runtime table
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown transportLaw type: " << name << nl
            << "Valid transportLaw types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // Construct and return the selected model
    return cstrIter()(name, dict, mesh, thermo);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

/*
 * read():
 * - Base no-op. Derived classes should:
 *     * extract required parameters from coeff()
 *     * validate ranges/units where appropriate
 */
void Foam::transportLaw::read()
{
    // No-op in base.
}


void Foam::transportLaw::calcMu
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& mu
) const
{
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        mu.cellField()[gpI] = calcMu(T.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        mu.faceField().minusValueAt(gpI) = calcMu(T.faceField().minusValue(gpI));
        mu.faceField().plusValueAt(gpI) = calcMu(T.faceField().plusValue(gpI));
    }
}


void Foam::transportLaw::calcMu
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& mu
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        mu[gpI] = calcMu(T[gpI]);
    }
}


void Foam::transportLaw::calcKappa
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& kappa
) const
{
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        kappa.cellField()[gpI] = calcKappa(T.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        kappa.faceField().minusValueAt(gpI) = calcKappa(T.faceField().minusValue(gpI));
        kappa.faceField().plusValueAt(gpI) = calcKappa(T.faceField().plusValue(gpI));
    }
}


void Foam::transportLaw::calcKappa
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& kappa
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        kappa[gpI] = calcKappa(T[gpI]);
    }
}


void Foam::transportLaw::calcPr
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& Pr
) const
{
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        Pr.cellField()[gpI] = calcPr(T.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        Pr.faceField().minusValueAt(gpI) = calcPr(T.faceField().minusValue(gpI));
        Pr.faceField().plusValueAt(gpI) = calcPr(T.faceField().plusValue(gpI));
    }
}


void Foam::transportLaw::calcPr
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& Pr
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        Pr[gpI] = calcPr(T[gpI]);
    }
}

// ************************************************************************* //

} // End namespace Foam
