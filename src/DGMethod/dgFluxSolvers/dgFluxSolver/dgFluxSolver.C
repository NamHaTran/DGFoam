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

#include "dgFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{

defineTypeNameAndDebug(dgFluxSolver, 0);
defineRunTimeSelectionTable(dgFluxSolver, dictionary);


Foam::dgFluxSolver::dgFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    name_(name),
    dict_(dict),
    mesh_(mesh)
{}


autoPtr<dgFluxSolver> Foam::dgFluxSolver::New
(
    const word& name,
    const word& fluxSolverType,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(fluxSolverType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgFluxSolver type: " << fluxSolverType << nl
            << "Valid types are: " << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dict, mesh);
}

void Foam::dgFluxSolver::makeONB
(
    const vector& n,
    vector& t1,
    vector& t2
) const
{
    // Construct t1 not colinear with n
    vector a = (mag(n.x()) < 0.9) ? vector(1,0,0) : vector(0,1,0);
    t1 = a - (a & n)*n;           // remove normal component
    scalar m = mag(t1) + SMALL;   // normalize with guard
    t1 /= m;
    t2 = n ^ t1;                  // right-handed basis
}

void Foam::dgFluxSolver::decomposeU
(
    const vector& U,
    const vector& n,
    vector& Un,
    vector& Ut1,
    vector& Ut2
) const
{
    // Build orthonormal basis (n, t1, t2)
    vector t1, t2;
    makeONB(n, t1, t2);

    // Components of U along basis directions
    const scalar Un_scalar  = U & n;     // projection on n
    const scalar Ut1_scalar = U & t1;    // projection on t1
    const scalar Ut2_scalar = U & t2;    // projection on t2

    // Convert back to vector form along the basis
    Un  = Un_scalar  * n;
    Ut1 = Ut1_scalar * t1;
    Ut2 = Ut2_scalar * t2;
}

} // End namespace Foam
// ************************************************************************* //

