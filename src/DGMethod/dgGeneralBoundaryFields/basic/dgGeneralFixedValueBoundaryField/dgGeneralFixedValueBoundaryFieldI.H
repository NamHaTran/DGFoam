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
    Foam::dgGeneralFixedValueBoundaryField

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralFixedValueBoundaryFieldI_H
#define dgGeneralFixedValueBoundaryFieldI_H

#include "dgGeneralBoundaryField.H"

namespace Foam
{

// * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

template<class Type>
dgGeneralFixedValueBoundaryField<Type>::dgGeneralFixedValueBoundaryField
(
    const word& name,
    const dgGeomMesh& dgMesh,
    const dictionary& dict
)
:
    dgGeneralBoundaryField<Type>(name, dgMesh, dict)
{
    dict.lookup("value") >> value_;
}

} // End namespace Foam

#endif
// ************************************************************************* //