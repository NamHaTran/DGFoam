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

#include "PerssonPeraire.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{

defineTypeNameAndDebug(PerssonPeraire, 0);
addToRunTimeSelectionTable(troubleCellDetector, PerssonPeraire, dictionary);


PerssonPeraire::PerssonPeraire
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    troubleCellDetector(dict, mesh)
{}


bool PerssonPeraire::detect(const label cellID) const
{
    const GaussField<scalar>& rhoGF = rhoGauss(cellID);
    const dofField<scalar>* rhoDofPtr = rhoGF.dofFieldPtr();

    if (!rhoDofPtr)
    {
        FatalErrorInFunction
            << "PerssonPeraire detector requires a density dofField."
            << abort(FatalError);
    }

    const List<scalar>& rhoDof = (*rhoDofPtr)[cellID].dof();

    if (rhoDof.size() <= 1)
    {
        return false;
    }

    const dgGeomCell& cell = *mesh_.cells()[cellID];
    const List<List<scalar>>& basis = cell.basis();
    const List<scalar>& weights = cell.weights();
    const List<scalar>& J3D = cell.J3D();

    const label highestMode = rhoDof.size() - 1;

    scalar numerator = Zero;
    scalar denominator = Zero;

    // Reconstruct the highest modal contribution and the full density field
    // at cell Gauss points, then compare their integrated energies.
    forAll(weights, gpI)
    {
        scalar rhoHigh = basis[gpI][highestMode]*rhoDof[highestMode];
        scalar rhoVal = Zero;

        forAll(rhoDof, dofI)
        {
            rhoVal += basis[gpI][dofI]*rhoDof[dofI];
        }

        const scalar w = weights[gpI]*J3D[gpI];

        numerator += w*sqr(rhoHigh);
        denominator += w*sqr(rhoVal);
    }

    if (denominator <= VSMALL)
    {
        return false;
    }

    // Preserve the threshold used by the legacy 2D implementation.
    const scalar sensor = numerator/denominator;
    const scalar threshold = 1.0/pow(scalar(mesh_.pOrder() + 1), 6.0);

    return sensor > threshold;
}

} // End namespace Foam

// ************************************************************************* //
