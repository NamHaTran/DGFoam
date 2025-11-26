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

Class
    Foam::dgInletOutletBoundaryField

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralInletOutletBoundaryFieldI_H
#define dgGeneralInletOutletBoundaryFieldI_H

#include "dgField.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

template<class Type>
dgGeneralInletOutletBoundaryField<Type>::dgGeneralInletOutletBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dictionary& dict
)
:
    dgGeneralBoundaryField<Type>(patch, dgMesh, dict),
    UName_(dict.lookupOrDefault("U", word("U"))),
    U_
    (
        dgMesh.getFvMesh().lookupObject<dgField<vector>>(UName_)
    )
{
    dict.readEntry("inletValue", inletValue_);
}

} // End namespace Foam

#endif
// ************************************************************************* //