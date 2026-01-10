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

#include "dgLaxFriedrichsFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dgLaxFriedrichsFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgLaxFriedrichsFluxSolver, dictionary);


// * * * * * * * * * * Constructor * * * * * * * * * * //

dgLaxFriedrichsFluxSolver::dgLaxFriedrichsFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgFluxSolver(name, dict, mesh),

    U_( mesh.getFvMesh().lookupObject<dgField<vector>>("U") ),
    a_( mesh.getFvMesh().lookupObject<dgField<scalar>>("a") ),
    scaleByMach_(false)
{
    read(dict);
}


// * * * * * * * * * * Read dict * * * * * * * * * * //

void dgLaxFriedrichsFluxSolver::read(const dictionary& dict)
{
    scaleByMach_ = dict.lookupOrDefault<bool>("scaleByMach", false);
}


// * * * * * * * * * * Dissipation coefficient * * * * * * * * * * //

scalar dgLaxFriedrichsFluxSolver::calcDissipationCoeff
(
    const vector& ULv,
    const vector& URv,
    const scalar aL,
    const scalar aR,
    const vector& n
) const
{
    const scalar UnL = (ULv & n);
    const scalar UnR = (URv & n);

    if (!scaleByMach_)
    {
        // Standard Rusanov/LF
        const scalar CL = mag(UnL) + aL;
        const scalar CR = mag(UnR) + aR;
        return max(CL, CR);
    }
    else
    {
        // Mach-scaled version
        scalar ML = (aL > SMALL ? mag(UnL)/aL : 0.0);
        scalar MR = (aR > SMALL ? mag(UnR)/aR : 0.0);

        const scalar invML = (ML > SMALL ? 1.0/ML : 1.0/SMALL);
        const scalar invMR = (MR > SMALL ? 1.0/MR : 1.0/SMALL);

        const scalar CL = mag(UnL) + aL*invML;
        const scalar CR = mag(UnR) + aR*invMR;

        return max(CL, CR);
    }
}


// * * * * * * * * * * Scalar U, vector F * * * * * * * * * * //

void dgLaxFriedrichsFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<vector>& F,
    const faceGaussField<scalar>& U
)
{
    const faceGaussField<vector>& UF = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& aF = a_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI=0; fI<nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI=0; gI<nGauss; ++gI)
        {
            const vector ULv = UF.minusValueOnFace(fI, gI);
            const vector URv = UF.plusValueOnFace(fI, gI);
            const scalar aL = aF.minusValueOnFace(fI, gI);
            const scalar aR = aF.plusValueOnFace(fI, gI);

            const scalar C = calcDissipationCoeff(ULv, URv, aL, aR, n);

            const vector FL = F.minusValueOnFace(fI, gI);
            const vector FR = F.plusValueOnFace(fI, gI);

            const scalar fL = (FL & n);
            const scalar fR = (FR & n);

            const scalar ULs = U.minusValueOnFace(fI, gI);
            const scalar URs = U.plusValueOnFace(fI, gI);

            const scalar fn = 0.5*(fL + fR) - 0.5*C*(URs - ULs);

            F.fluxOnFace(fI, gI) = fn*n;
        }
    }
}


// * * * * * * * * * * Vector U, tensor F * * * * * * * * * * //

void dgLaxFriedrichsFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>& U
)
{
    const faceGaussField<vector>& UF = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& aF = a_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI=0; fI<nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI=0; gI<nGauss; ++gI)
        {
            const vector ULv = UF.minusValueOnFace(fI, gI);
            const vector URv = UF.plusValueOnFace(fI, gI);
            const scalar aL = aF.minusValueOnFace(fI, gI);
            const scalar aR = aF.plusValueOnFace(fI, gI);

            const scalar C = calcDissipationCoeff(ULv, URv, aL, aR, n);

            const tensor FL = F.minusValueOnFace(fI, gI);
            const tensor FR = F.plusValueOnFace(fI, gI);

            const vector fL = (FL & n);
            const vector fR = (FR & n);

            const vector ULs = U.minusValueOnFace(fI, gI);
            const vector URs = U.plusValueOnFace(fI, gI);

            const vector fn = 0.5*(fL + fR) - 0.5*C*(URs - ULs);

            F.fluxOnFace(fI, gI) = fn;   // already vector flux
        }
    }
}

} // End namespace Foam

// ************************************************************************* //

