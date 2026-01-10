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

#ifndef Foam_dgInterfaceBoundaryFieldI_H
#define Foam_dgInterfaceBoundaryFieldI_H

// * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::dgInterfaceBoundaryField<Type>::dgInterfaceBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh
)
:
    patch_(patch),
    dgMesh_(dgMesh)
{}


// * * * * * * * * * * * * * Factory Method  * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::dgInterfaceBoundaryField<Type>>
Foam::dgInterfaceBoundaryField<Type>::New
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh
)
{
    const word bcType(patch.type());
    
    auto cstrIter = fvPatchConstructorTablePtr_->find(bcType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown dgInterfaceBoundaryField type: " << bcType << nl
            << "Valid dgInterfaceBoundaryField types are: "
            << dgInterfaceBoundaryField<Type>::fvPatchConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // Call constructor supplied in table
    return cstrIter()(patch, dgMesh);
}

#endif
// ************************************************************************* //