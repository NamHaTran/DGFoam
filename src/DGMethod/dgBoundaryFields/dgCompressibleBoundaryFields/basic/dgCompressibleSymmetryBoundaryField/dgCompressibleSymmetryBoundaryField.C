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

#include "dgCompressibleSymmetryBoundaryField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleSymmetryBoundaryField, 0);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressibleSymmetryBoundaryField,
    dictionary
);

dgCompressibleSymmetryBoundaryField::dgCompressibleSymmetryBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict)
{}


void dgCompressibleSymmetryBoundaryField::updateGhostState
(
    const label,
    const label,
    const label,
    const vector& n,
    const scalar rhoMinus,
    const vector& rhoUMinus,
    const scalar EMinus,
    scalar& rhoPlus,
    vector& rhoUPlus,
    scalar& EPlus
) const
{
    rhoPlus  = rhoMinus;
    rhoUPlus = rhoUMinus - 2.0*(n & rhoUMinus)*n;
    EPlus    = EMinus;
}


void dgCompressibleSymmetryBoundaryField::updateBCValue
(
    const label,
    const label,
    const label,
    const vector& n,
    const scalar rhoMinus,
    const vector& rhoUMinus,
    const scalar EMinus,
    scalar& rhoBC,
    vector& rhoUBC,
    scalar& EBC
) const
{
    rhoBC  = rhoMinus;
    rhoUBC = rhoUMinus - (n & rhoUMinus)*n;
    EBC    = EMinus;
}


void dgCompressibleSymmetryBoundaryField::checkPatchType() const
{
    if (this->patch_.type() != "symmetry")
    {
        FatalErrorInFunction
            << "Boundary condition " << this->type() << " can only "
            << "be applied to patch type:\n"
            << "    symmetry\n"
            << "but patch " << this->patch_.name()
            << " is of type " << this->patch_.type() << nl
            << exit(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
