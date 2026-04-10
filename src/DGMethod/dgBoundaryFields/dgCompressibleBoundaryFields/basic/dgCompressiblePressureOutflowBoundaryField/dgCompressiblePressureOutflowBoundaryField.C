/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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

#include "dgCompressiblePressureOutflowBoundaryField.H"
#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressiblePressureOutflowBoundaryField, 0);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressiblePressureOutflowBoundaryField,
    dictionary
);

dgCompressiblePressureOutflowBoundaryField::
dgCompressiblePressureOutflowBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict),
    pValue_(Zero),
    subsonicMachThreshold_(0.99)
{
    const dictionary& coeffDict =
        dict.subDict("compressiblePressureOutflowCoeff");

    pValue_ = max(coeffDict.get<scalar>("pValue"), scalar(SMALL));
    subsonicMachThreshold_ =
        max
        (
            coeffDict.lookupOrDefault<scalar>("subsonicMachThreshold", 0.99),
            scalar(0)
        );
}


void dgCompressiblePressureOutflowBoundaryField::updateGhostState
(
    const label cellID,
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
    const scalar rhoMinusSafe = max(rhoMinus, scalar(SMALL));
    const vector UMinus = rhoUMinus/rhoMinusSafe;
    const scalar kMinus = 0.5*magSqr(UMinus);
    const scalar heMinus = EMinus/rhoMinusSafe - kMinus;
    const scalar aMinus =
        max(thermo_.calcSpeedOfSoundFromRhoHe(cellID, rhoMinusSafe, heMinus), scalar(SMALL));
    const scalar machMinus = mag(n & UMinus)/aMinus;

    if (machMinus >= subsonicMachThreshold_)
    {
        rhoPlus = rhoMinus;
        rhoUPlus = rhoUMinus;
        EPlus = EMinus;
        return;
    }

    rhoPlus = rhoMinus;
    rhoUPlus = rhoUMinus;
    EPlus = pressureVelocityToConservativeEnergy(pValue_, rhoMinusSafe, UMinus);
}


void dgCompressiblePressureOutflowBoundaryField::checkPatchType() const
{
    if (this->patch_.type() != "patch")
    {
        FatalErrorInFunction
            << "Boundary condition " << this->type() << " can only "
            << "be applied to patch type:\n"
            << "    patch\n"
            << "but patch " << this->patch_.name()
            << " is of type " << this->patch_.type() << nl
            << exit(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
