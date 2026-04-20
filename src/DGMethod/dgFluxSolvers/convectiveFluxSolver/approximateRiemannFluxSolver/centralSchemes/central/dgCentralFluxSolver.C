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

#include "dgCentralFluxSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dgCentralFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgCentralFluxSolver, dictionary);


// * * * * * * * * * Scalar U, vector F  * * * * * * * * * //

void dgCentralFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<vector>& F,
    const faceGaussField<scalar>& U
)
{
    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI=0; fI<nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI=0; gI<nGauss; ++gI)
        {
            // Physical fluxes (vector)
            const vector FL = F.minusValueOnFace(fI, gI);
            const vector FR = F.plusValueOnFace(fI, gI);

            // Project to normal direction → scalar
            const scalar fL = (FL & n);
            const scalar fR = (FR & n);

            // Central scalar flux
            const scalar fC = 0.5*(fL + fR);

            // Back to Cartesian: vector = scalar * n
            F.fluxOnFace(fI, gI) = fC*n;
        }
    }
}


// * * * * * * * * * Vector U, tensor F  * * * * * * * * * //

void dgCentralFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>& U
)
{
    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI=0; fI<nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI=0; gI<nGauss; ++gI)
        {
            // Physical fluxes (tensor)
            const tensor FL = F.minusValueOnFace(fI, gI);
            const tensor FR = F.plusValueOnFace(fI, gI);

            // Project tensor → vector by dot normal
            const vector fL = (FL & n);
            const vector fR = (FR & n);

            // Central vector flux
            const vector fC = 0.5*(fL + fR);

            // Final flux is already a vector (no projection needed)
            F.fluxOnFace(fI, gI) = fC;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //

