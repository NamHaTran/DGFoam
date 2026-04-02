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

#include "energy.H"
#include "addToRunTimeSelectionTable.H"
#include "dgExpr.H"

namespace Foam
{

defineTypeNameAndDebug(energy, 0);
defineRunTimeSelectionTable(energy, dictionary);

word energy::heName(const heType heTypeValue)
{
    switch (heTypeValue)
    {
        case heType::internalEnergy:
            return "e";

        case heType::enthalpy:
            return "h";
    }

    FatalErrorInFunction
        << "Unknown heType value." << exit(FatalError);

    return word::null;
}


energy::energy
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo,
    const heType heTypeValue
)
:
    name_(name),
    coeffs_(dict),
    mesh_(mesh),
    thermo_(thermo),
    he_(energy::heName(heTypeValue)),
    heType_(heTypeValue)
{}

autoPtr<energy> energy::New(const word& name, const dictionary& dict, const dgGeomMesh& mesh, const thermoLaw& thermo)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown energy model type: " << name << nl
            << "Valid types are: " << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dict, mesh, thermo);
}


void energy::calcEnthalpy
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
                return calcEnthalpy(TValue);
            },
            dg::expr(T)
        )
    );
}


void energy::calcEnthalpy
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& h
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        h[gpI] = calcEnthalpy(T[gpI]);
    }
}


void energy::calcHe
(
    const label,
    const GaussField<scalar>& T,
    GaussField<scalar>& he
) const
{
    dg::assign
    (
        he,
        dg::map
        (
            [this](const scalar TValue)
            {
                return calcHe(TValue);
            },
            dg::expr(T)
        )
    );
}


void energy::calcHe
(
    const boundaryGaussField<scalar>& T,
    boundaryGaussField<scalar>& he
) const
{
    for (label gpI = 0; gpI < T.size(); ++gpI)
    {
        he[gpI] = calcHe(T[gpI]);
    }
}


void energy::calcTfromHe
(
    const label,
    const GaussField<scalar>& he,
    GaussField<scalar>& T
) const
{
    dg::assign
    (
        T,
        dg::map
        (
            [this](const scalar heValue)
            {
                return calcTfromHe(heValue);
            },
            dg::expr(he)
        )
    );
}


void energy::calcTfromHe
(
    const boundaryGaussField<scalar>& he,
    boundaryGaussField<scalar>& T
) const
{
    for (label gpI = 0; gpI < he.size(); ++gpI)
    {
        T[gpI] = calcTfromHe(he[gpI]);
    }
}

} // namespace Foam

// ************************************************************************* //
