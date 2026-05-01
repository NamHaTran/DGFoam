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

#ifndef Foam_EulerFieldStateI_H
#define Foam_EulerFieldStateI_H

#include "EulerFieldState.H"
#include "error.H"
#include "IOstreams.H"
#include "pTraits.H"

namespace Foam
{

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

template<class Type>
EulerFieldState<Type>::EulerFieldState
(
    dgField<Type>& field,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    U_(field),
    UOld_
    (
        IOobject
        (
            field.name() + "Old",
            mesh.getFvMesh().time().timeName(),
            mesh.getFvMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        field.hasDof()
    ),
    R_()
{
    // Allocate buffers once, then start from a clean saved state.
    allocateResidual();
    saveOldField();
}


// * * * * * * * * * * * * * Private Helpers * * * * * * * * * * * * * * * * //

template<class Type>
void EulerFieldState<Type>::allocateResidual()
{
    // Residual storage follows the mesh DoF layout cell by cell.
    R_.setSize(mesh_.nCells());

    forAll(R_, cellI)
    {
        const label nDof = mesh_.cells()[cellI]->nDof();
        R_[cellI].setSize(nDof, pTraits<Type>::zero);
    }
}


template<class Type>
void EulerFieldState<Type>::copyField
(
    dgField<Type>& dst,
    const dgField<Type>& src
)
{
    // Copy the full DG field state into an internal buffer field.
    if (dst.hasDof() != src.hasDof())
    {
        FatalErrorInFunction
            << "Incompatible dgField copy between " << src.name()
            << " and " << dst.name() << ": hasDof mismatch."
            << abort(FatalError);
    }

    if (dst.nCells() != src.nCells())
    {
        FatalErrorInFunction
            << "Incompatible dgField copy between " << src.name()
            << " and " << dst.name() << ": cell-count mismatch."
            << abort(FatalError);
    }

    if (src.hasDof())
    {
        forAll(dst.gaussFields(), cellI)
        {
            dst.dof()[cellI] = src.dof()[cellI];
            dst.gaussFields()[cellI] = src.gaussFields()[cellI];
            dst.dof().updateCellDof(cellI);
        }
    }
    else
    {
        forAll(dst.gaussFields(), cellI)
        {
            dst.gaussFields()[cellI] = src.gaussFields()[cellI];
        }
    }
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class Type>
void EulerFieldState<Type>::saveOldField()
{
    // Freeze U^n before the scheme starts updating U.
    copyField(UOld_, U_);
}


template<class Type>
void EulerFieldState<Type>::restoreOldField()
{
    // Rebuild the active field from the saved step start.
    copyField(U_, UOld_);
}


template<class Type>
void EulerFieldState<Type>::clearResidual()
{
    // Each new step starts from zero residual.
    forAll(R_, cellI)
    {
        forAll(R_[cellI], dofI)
        {
            R_[cellI][dofI] = pTraits<Type>::zero;
        }
    }
}


template<class Type>
void EulerFieldState<Type>::reset()
{
    // Residual storage is reused and cleared just before assembly.
    restoreOldField();
}

} // End namespace Foam

#endif

// ************************************************************************* //
