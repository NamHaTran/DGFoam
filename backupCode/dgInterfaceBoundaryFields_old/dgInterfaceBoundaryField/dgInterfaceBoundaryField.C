/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DGFoam: Discontinuous Galerkin CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | GPU-friendly CFD solver framework
     \\/     M anipulation  |
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgInterfaceBoundaryField.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members  * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(dgInterfaceBoundaryField<scalar>, 0);
defineNamedTemplateTypeNameAndDebug(dgInterfaceBoundaryField<vector>, 0);
defineNamedTemplateTypeNameAndDebug(dgInterfaceBoundaryField<tensor>, 0);

defineTemplateRunTimeSelectionTable(dgInterfaceBoundaryField<scalar>, fvPatch);
defineTemplateRunTimeSelectionTable(dgInterfaceBoundaryField<vector>, fvPatch);
defineTemplateRunTimeSelectionTable(dgInterfaceBoundaryField<tensor>, fvPatch);

} // End namespace Foam

// ************************************************************************* //

