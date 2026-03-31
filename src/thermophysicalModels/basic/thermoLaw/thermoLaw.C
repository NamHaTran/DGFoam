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
    for (label gpI = 0; gpI < e.cellField().size(); ++gpI)
    {
        T.cellField()[gpI] = calcTFromInternalE(e.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < e.faceField().nGauss(); ++gpI)
    {
        T.faceField().minusValueAt(gpI) =
            calcTFromInternalE(e.faceField().minusValue(gpI));
        T.faceField().plusValueAt(gpI) =
            calcTFromInternalE(e.faceField().plusValue(gpI));
    }
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
    for (label gpI = 0; gpI < h.cellField().size(); ++gpI)
    {
        T.cellField()[gpI] = calcTFromH(h.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < h.faceField().nGauss(); ++gpI)
    {
        T.faceField().minusValueAt(gpI) =
            calcTFromH(h.faceField().minusValue(gpI));
        T.faceField().plusValueAt(gpI) =
            calcTFromH(h.faceField().plusValue(gpI));
    }
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
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        Cp.cellField()[gpI] = calcCp(T.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        Cp.faceField().minusValueAt(gpI) = calcCp(T.faceField().minusValue(gpI));
        Cp.faceField().plusValueAt(gpI) = calcCp(T.faceField().plusValue(gpI));
    }
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
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        Cv.cellField()[gpI] = calcCv(T.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        Cv.faceField().minusValueAt(gpI) = calcCv(T.faceField().minusValue(gpI));
        Cv.faceField().plusValueAt(gpI) = calcCv(T.faceField().plusValue(gpI));
    }
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
    for (label gpI = 0; gpI < Cp.cellField().size(); ++gpI)
    {
        gamma.cellField()[gpI] =
            calcGamma(Cp.cellField()[gpI], Cv.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < Cp.faceField().nGauss(); ++gpI)
    {
        gamma.faceField().minusValueAt(gpI) = calcGamma
        (
            Cp.faceField().minusValue(gpI),
            Cv.faceField().minusValue(gpI)
        );

        gamma.faceField().plusValueAt(gpI) = calcGamma
        (
            Cp.faceField().plusValue(gpI),
            Cv.faceField().plusValue(gpI)
        );
    }
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
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        e.cellField()[gpI] = calcInternalE(T.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        e.faceField().minusValueAt(gpI) =
            calcInternalE(T.faceField().minusValue(gpI));
        e.faceField().plusValueAt(gpI) =
            calcInternalE(T.faceField().plusValue(gpI));
    }
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
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        h.cellField()[gpI] = calcH(T.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        h.faceField().minusValueAt(gpI) = calcH(T.faceField().minusValue(gpI));
        h.faceField().plusValueAt(gpI) = calcH(T.faceField().plusValue(gpI));
    }
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
    for (label gpI = 0; gpI < T.cellField().size(); ++gpI)
    {
        a.cellField()[gpI] =
            calcSpeedOfSound(T.cellField()[gpI], gamma.cellField()[gpI]);
    }

    for (label gpI = 0; gpI < T.faceField().nGauss(); ++gpI)
    {
        a.faceField().minusValueAt(gpI) = calcSpeedOfSound
        (
            T.faceField().minusValue(gpI),
            gamma.faceField().minusValue(gpI)
        );

        a.faceField().plusValueAt(gpI) = calcSpeedOfSound
        (
            T.faceField().plusValue(gpI),
            gamma.faceField().plusValue(gpI)
        );
    }
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
