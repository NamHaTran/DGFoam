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

#include "dgGeneralBoundaryField.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members  * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(dgGeneralBoundaryField<scalar>, 0);
defineNamedTemplateTypeNameAndDebug(dgGeneralBoundaryField<vector>, 0);
defineNamedTemplateTypeNameAndDebug(dgGeneralBoundaryField<tensor>, 0);

defineTemplateRunTimeSelectionTable(dgGeneralBoundaryField<scalar>, dictionary);
defineTemplateRunTimeSelectionTable(dgGeneralBoundaryField<vector>, dictionary);
defineTemplateRunTimeSelectionTable(dgGeneralBoundaryField<tensor>, dictionary);

} // End namespace Foam

// ************************************************************************* //

