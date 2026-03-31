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

#include "dgCompressibleFixedValueBoundaryField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dgCompressibleFixedValueBoundaryField, 0);
addToRunTimeSelectionTable
(
    dgCompressibleBoundaryField,
    dgCompressibleFixedValueBoundaryField,
    dictionary
);

dgCompressibleFixedValueBoundaryField::dgCompressibleFixedValueBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dgThermoConservative& thermo,
    const dictionary& dict
)
:
    dgCompressibleBoundaryField(patch, dgMesh, thermo, dict),
    rhoValue_(Zero),
    rhoUValue_(Zero),
    EValue_(Zero)
{
    const dictionary& coeffDict = dict.subDict("compressibleFixedValueCoeff");

    const scalar pValue = coeffDict.get<scalar>("pValue");
    const scalar TValue = coeffDict.get<scalar>("TValue");
    const vector UValue = coeffDict.get<vector>("UValue");

    primitiveToConservative(pValue, TValue, UValue, rhoValue_, rhoUValue_, EValue_);
}


void dgCompressibleFixedValueBoundaryField::updateGhostState
(
    const label,
    const label,
    const label,
    const vector&,
    const scalar,
    const vector&,
    const scalar,
    scalar& rhoPlus,
    vector& rhoUPlus,
    scalar& EPlus
) const
{
    rhoPlus  = rhoValue_;
    rhoUPlus = rhoUValue_;
    EPlus    = EValue_;
}


void dgCompressibleFixedValueBoundaryField::checkPatchType() const
{
    if
    (
        this->patch_.type() != "patch"
     && this->patch_.type() != "wall"
    )
    {
        FatalErrorInFunction
            << "Boundary condition " << this->type() << " can only "
            << "be applied to patch types:\n"
            << "    patch\n"
            << "    wall\n"
            << "but patch " << this->patch_.name()
            << " is of type " << this->patch_.type() << nl
            << exit(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
