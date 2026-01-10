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

#ifndef Foam_dgProcessorBoundaryManagerI_H
#define Foam_dgProcessorBoundaryManagerI_H

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::dgProcessorBoundaryManager<Type>::dgProcessorBoundaryManager
(
    const dgGeomMesh& mesh
)
:
    mesh_(mesh)
{
    const fvBoundaryMesh& patches = mesh_.getFvMesh().boundary();
    const label nPatches = patches.size();

    for (label patchI = 0; patchI < nPatches; ++patchI)
    {
        const fvPatch& patch = patches[patchI];

        if (dgBoundaryHelper::isInterfacePatch(patch))
        {
            bConditions_.append
            (
                new dgProcessorBoundaryField<Type>(patch, mesh_)
            );
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::dgProcessorBoundaryManager<Type>::update
(
    dgField<Type>& field
)
{
    if (Pstream::parRun())
    {
        // pack and send
        forAll(bConditions_, patchI)
        {
            bConditions_[patchI].sendAll(field);
        }

        // Wait for all sends to complete
        UPstream::waitRequests(-1);

        // recv and unpack
        forAll(bConditions_, patchI)
        {
            bConditions_[patchI].recvAll(field);
        }
    }
}

#endif
// ************************************************************************* //

