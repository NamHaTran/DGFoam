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

#ifndef Foam_SSPRK4FieldStateI_H
#define Foam_SSPRK4FieldStateI_H

#include "SSPRK4FieldState.H"
#include "error.H"
#include "pTraits.H"

namespace Foam
{

template<class Type>
SSPRK4FieldState<Type>::SSPRK4FieldState
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
    UStages_(),
    RStages_()
{
    // SSPRK4 needs one old field, four stage fields, and five residuals.
    allocateStageFields();
    allocateResiduals();
    saveOldField();
    clearResiduals();
}


template<class Type>
void SSPRK4FieldState<Type>::allocateStageFields()
{
    // Allocate explicit storage for all intermediate SSP states.
    UStages_.setSize(nStoredStages_);

    for (label stageFieldI = 0; stageFieldI < nStoredStages_; ++stageFieldI)
    {
        UStages_.set
        (
            stageFieldI,
            new dgField<Type>
            (
                IOobject
                (
                    U_.name() + "Stage" + Foam::name(stageFieldI + 1),
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                U_.hasDof()
            )
        );
    }
}


template<class Type>
void SSPRK4FieldState<Type>::allocateResiduals()
{
    // Residual buffers follow the mesh DoF layout for each stage.
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
void SSPRK4FieldState<Type>::copyField
(
    dgField<Type>& dst,
    const dgField<Type>& src
)
{
    // Copy whole-field DG data into another buffer field.
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


template<class Type>
void SSPRK4FieldState<Type>::saveOldField()
{
    // Save U^n once before entering the SSPRK4 loop.
    copyField(UOld_, U_);
}


template<class Type>
void SSPRK4FieldState<Type>::restoreOldField()
{
    // Restore the active field from the saved step start.
    copyField(U_, UOld_);
}


template<class Type>
void SSPRK4FieldState<Type>::saveStageField(const label stageFieldI)
{
    // Store one intermediate SSP stage field.
    copyField(UStages_[stageFieldI], U_);
}


template<class Type>
void SSPRK4FieldState<Type>::restoreStageField(const label stageFieldI)
{
    // Put one stored stage back into the active field.
    copyField(U_, UStages_[stageFieldI]);
}


template<class Type>
void SSPRK4FieldState<Type>::clearResiduals()
{
    // Clear all residual buffers together.
    forAll(RStages_, stageI)
    {
        clearResidual(stageI);
    }
}


template<class Type>
void SSPRK4FieldState<Type>::clearResidual(const label stageI)
{
    if (stageI < 0 || stageI >= nStages_)
    {
        FatalErrorInFunction
            << "Invalid SSPRK4 stage index " << stageI
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
void SSPRK4FieldState<Type>::reset()
{
    // Reset restores U and refreshes all stored stage fields.
    restoreOldField();
    for (label stageFieldI = 0; stageFieldI < nStoredStages_; ++stageFieldI)
    {
        copyField(UStages_[stageFieldI], UOld_);
    }
    clearResiduals();
}

} // End namespace Foam

#endif

// ************************************************************************* //
