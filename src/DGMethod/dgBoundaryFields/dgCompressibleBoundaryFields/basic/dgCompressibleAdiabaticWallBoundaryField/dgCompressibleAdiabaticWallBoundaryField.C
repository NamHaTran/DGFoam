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

#include "dgCompressibleAdiabaticWallBoundaryField.H"
#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleAdiabaticWallBoundaryField, 0);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressibleAdiabaticWallBoundaryField,
    dictionary
);

dgCompressibleAdiabaticWallBoundaryField::
dgCompressibleAdiabaticWallBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict)
{}


void dgCompressibleAdiabaticWallBoundaryField::updateGhostState
(
    const label,
    const label,
    const label,
    const vector&,
    const scalar rhoMinus,
    const vector& rhoUMinus,
    const scalar EMinus,
    scalar& rhoPlus,
    vector& rhoUPlus,
    scalar& EPlus
) const
{
    const vector UMinus = rhoUMinus/rhoMinus;
    const scalar kMinus = 0.5*magSqr(UMinus);
    const scalar heMinus = EMinus/rhoMinus - kMinus;
    const scalar TMinus = thermo_.calcTemperatureFromRhoHe(rhoMinus, heMinus);
    const scalar pMinus = thermo_.calcPressureFromRhoHe(rhoMinus, heMinus);

    primitiveToConservative
    (
        pMinus,
        TMinus,
        vector(Zero),
        rhoPlus,
        rhoUPlus,
        EPlus
    );
}


void dgCompressibleAdiabaticWallBoundaryField::checkPatchType() const
{
    if (this->patch_.type() != "wall")
    {
        FatalErrorInFunction
            << "Boundary condition " << this->type() << " can only "
            << "be applied to patch type:\n"
            << "    wall\n"
            << "but patch " << this->patch_.name()
            << " is of type " << this->patch_.type() << nl
            << exit(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
