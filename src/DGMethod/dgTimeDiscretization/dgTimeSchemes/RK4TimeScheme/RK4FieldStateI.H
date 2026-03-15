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

#ifndef Foam_RK4FieldStateI_H
#define Foam_RK4FieldStateI_H

#include "RK4FieldState.H"
#include "error.H"
#include "pTraits.H"

namespace Foam
{

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

template<class Type>
RK4FieldState<Type>::RK4FieldState
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
    RStages_()
{
    // RK4 needs one saved field and four residual slots.
    allocateResiduals();
    saveOldField();
    clearResiduals();
}


// * * * * * * * * * * * * * Private Helpers * * * * * * * * * * * * * * * * //

template<class Type>
void RK4FieldState<Type>::allocateResiduals()
{
    // Each stage keeps its own residual on the mesh DoF layout.
    RStages_.setSize(nStages_);

    forAll(RStages_, stageI)
    {
        RStages_[stageI].setSize(mesh_.nCells());

        forAll(RStages_[stageI], cellI)
        {
            const label nDof = mesh_.cells()[cellI]->nDof();
            RStages_[stageI][cellI].setSize(nDof, pTraits<Type>::zero);
        }
    }
}


template<class Type>
void RK4FieldState<Type>::copyField
(
    dgField<Type>& dst,
    const dgField<Type>& src
)
{
    // Copy whole-field data into another DG buffer field.
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
        }

        dst.dof().updateFoamFields();
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
void RK4FieldState<Type>::saveOldField()
{
    // Save U^n once before stage updates start.
    copyField(UOld_, U_);
}


template<class Type>
void RK4FieldState<Type>::restoreOldField()
{
    // Restore the working field from the saved step start.
    copyField(U_, UOld_);
}


template<class Type>
void RK4FieldState<Type>::clearResiduals()
{
    // Clear all stage residuals together.
    forAll(RStages_, stageI)
    {
        clearResidual(stageI);
    }
}


template<class Type>
void RK4FieldState<Type>::clearResidual(const label stageI)
{
    // Clear only the stage that is about to be assembled.
    if (stageI < 0 || stageI >= nStages_)
    {
        FatalErrorInFunction
            << "Invalid RK4 stage index " << stageI
            << ". Valid stage indices are 0.." << (nStages_ - 1) << "."
            << abort(FatalError);
    }

    forAll(RStages_[stageI], cellI)
    {
        forAll(RStages_[stageI][cellI], dofI)
        {
            RStages_[stageI][cellI][dofI] = pTraits<Type>::zero;
        }
    }
}


template<class Type>
void RK4FieldState<Type>::reset()
{
    // Reset both the active field and all stage storage.
    restoreOldField();
    clearResiduals();
}

} // End namespace Foam

#endif

// ************************************************************************* //
