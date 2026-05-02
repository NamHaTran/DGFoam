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

#include "dgSSPRK3TimeScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "pTraits.H"

namespace Foam
{

defineTypeNameAndDebug(dgSSPRK3TimeScheme, 0);
addToRunTimeSelectionTable(dgTimeScheme, dgSSPRK3TimeScheme, dictionary);

template<class Type>
void copyCell(dgField<Type>& dst, const dgField<Type>& src, const label cellID)
{
    if (dst.hasDof() != src.hasDof())
    {
        FatalErrorInFunction
            << "Incompatible dgField cell copy between " << src.name()
            << " and " << dst.name() << ": hasDof mismatch."
            << abort(FatalError);
    }

    dst.gaussFields()[cellID] = src.gaussFields()[cellID];

    if (src.hasDof())
    {
        dst.dof()[cellID] = src.dof()[cellID];
    }
}


template<class Type>
void saveOldCell(SSPRK3FieldState<Type>& state, const label cellID)
{
    copyCell(state.UOld(), state.U(), cellID);
}


template<class Type>
void restoreOldCell(SSPRK3FieldState<Type>& state, const label cellID)
{
    copyCell(state.U(), state.UOld(), cellID);
}


template<class Type>
void clearResidualCell
(
    SSPRK3FieldState<Type>& state,
    const label stageI,
    const label cellID
)
{
    List<Type>& RCell = state.R(stageI, cellID);

    forAll(RCell, dofI)
    {
        RCell[dofI] = pTraits<Type>::zero;
    }
}


template<class Type>
void syncUpdatedCell(dgField<Type>& U, const label cellID)
{
    U.dof().updateCellDof(cellID);
}


template<class Type>
void updateSSPRK3Stage
(
    SSPRK3FieldState<Type>& state,
    const dgGeomMesh& mesh,
    const label stageI,
    const label cellID,
    const scalar dt
)
{
    dgField<Type>& U = state.U();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot update SSPRK3 state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    List<Type>& UCellDof = U.dof()[cellID].dof();
    const List<Type>& UOldCellDof = state.UOld().dof()[cellID].dof();
    const List<scalar>& massDiag = mesh.cells()[cellID]->massMatrixDiag();

    switch (stageI)
    {
        case 0:
        {
            const List<Type>& R0 = state.R(0, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    UOldCellDof[dofI] + dt*(R0[dofI]/massDiag[dofI]);
            }

            syncUpdatedCell(U, cellID);
            break;
        }

        case 1:
        {
            const List<Type>& R1 = state.R(1, cellID);

            forAll(UCellDof, dofI)
            {
                const Type U1CellDof = UCellDof[dofI];

                UCellDof[dofI] =
                    0.75*UOldCellDof[dofI]
                  + 0.25*(U1CellDof + dt*(R1[dofI]/massDiag[dofI]));
            }

            syncUpdatedCell(U, cellID);
            break;
        }

        case 2:
        {
            const List<Type>& R2 = state.R(2, cellID);

            forAll(UCellDof, dofI)
            {
                const Type U2CellDof = UCellDof[dofI];

                UCellDof[dofI] =
                    (1.0/3.0)*UOldCellDof[dofI]
                  + (2.0/3.0)*(U2CellDof + dt*(R2[dofI]/massDiag[dofI]));
            }

            syncUpdatedCell(U, cellID);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Invalid SSPRK3 stage index " << stageI
                << ". Valid stage indices are 0..2."
                << abort(FatalError);
        }
    }
}


Foam::dgSSPRK3TimeScheme::dgSSPRK3TimeScheme
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgTimeScheme(name, dict, mesh),
    scalarNames_(),
    scalarStates_(),
    vectorNames_(),
    vectorStates_()
{
    read(dict);
}


void Foam::dgSSPRK3TimeScheme::read(const dictionary& dict)
{
    dict_ = dict;
}


label Foam::dgSSPRK3TimeScheme::findScalarIndex(const word& fieldName) const
{
    forAll(scalarNames_, i)
    {
        if (scalarNames_[i] == fieldName)
        {
            return i;
        }
    }

    return -1;
}


label Foam::dgSSPRK3TimeScheme::findVectorIndex(const word& fieldName) const
{
    forAll(vectorNames_, i)
    {
        if (vectorNames_[i] == fieldName)
        {
            return i;
        }
    }

    return -1;
}


void Foam::dgSSPRK3TimeScheme::validateStageIndex(const label stageI) const
{
    if (stageI < 0 || stageI >= nStages())
    {
        FatalErrorInFunction
            << "SSPRK3 scheme supports only stage indices 0.." << (nStages() - 1)
            << ", but got stageI = " << stageI << "."
            << abort(FatalError);
    }
}


void Foam::dgSSPRK3TimeScheme::validateStageIndex() const
{
    validateStageIndex(stage());
}


void Foam::dgSSPRK3TimeScheme::registerField(dgField<scalar>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register scalar field " << U.name()
            << " to dgSSPRK3TimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findScalarIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Scalar field " << U.name()
            << " is already registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    const label newI = scalarNames_.size();
    scalarNames_.setSize(newI + 1);
    scalarNames_[newI] = U.name();
    scalarStates_.setSize(newI + 1);
    scalarStates_.set(newI, new SSPRK3FieldState<scalar>(U, mesh_));
}


void Foam::dgSSPRK3TimeScheme::registerField(dgField<vector>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register vector field " << U.name()
            << " to dgSSPRK3TimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findVectorIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Vector field " << U.name()
            << " is already registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    const label newI = vectorNames_.size();
    vectorNames_.setSize(newI + 1);
    vectorNames_[newI] = U.name();
    vectorStates_.setSize(newI + 1);
    vectorStates_.set(newI, new SSPRK3FieldState<vector>(U, mesh_));
}


label Foam::dgSSPRK3TimeScheme::nStages() const
{
    return 3;
}


dgField<scalar>& Foam::dgSSPRK3TimeScheme::scalarField(const word& fieldName)
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


const dgField<scalar>& Foam::dgSSPRK3TimeScheme::scalarField(const word& fieldName) const
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


dgField<vector>& Foam::dgSSPRK3TimeScheme::vectorField(const word& fieldName)
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


const dgField<vector>& Foam::dgSSPRK3TimeScheme::vectorField(const word& fieldName) const
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


List<List<scalar>>& Foam::dgSSPRK3TimeScheme::scalarResidual
(
    const word& fieldName,
    const label stageI
)
{
    validateStageIndex(stageI);

    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R(stageI);
}


const List<List<scalar>>& Foam::dgSSPRK3TimeScheme::scalarResidual
(
    const word& fieldName,
    const label stageI
) const
{
    validateStageIndex(stageI);

    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R(stageI);
}


List<List<vector>>& Foam::dgSSPRK3TimeScheme::vectorResidual
(
    const word& fieldName,
    const label stageI
)
{
    validateStageIndex(stageI);

    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R(stageI);
}


const List<List<vector>>& Foam::dgSSPRK3TimeScheme::vectorResidual
(
    const word& fieldName,
    const label stageI
) const
{
    validateStageIndex(stageI);

    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgSSPRK3TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R(stageI);
}


void Foam::dgSSPRK3TimeScheme::saveOldFields(const label cellID)
{
    forAll(scalarStates_, i)
    {
        saveOldCell(scalarStates_[i], cellID);
    }

    forAll(vectorStates_, i)
    {
        saveOldCell(vectorStates_[i], cellID);
    }
}


void Foam::dgSSPRK3TimeScheme::restoreOldFields(const label cellID)
{
    forAll(scalarStates_, i)
    {
        restoreOldCell(scalarStates_[i], cellID);
    }

    forAll(vectorStates_, i)
    {
        restoreOldCell(vectorStates_[i], cellID);
    }
}


void Foam::dgSSPRK3TimeScheme::saveStageFields(const label cellID)
{
    // SSPRK3 keeps only U^n; the current U already carries intermediate stages.
}


void Foam::dgSSPRK3TimeScheme::restoreStageFields(const label cellID)
{
    // SSPRK3 does not persist intermediate stage fields.
}


void Foam::dgSSPRK3TimeScheme::clearResiduals(const label cellID)
{
    validateStageIndex();

    forAll(scalarStates_, i)
    {
        clearResidualCell(scalarStates_[i], stage(), cellID);
    }

    forAll(vectorStates_, i)
    {
        clearResidualCell(vectorStates_[i], stage(), cellID);
    }
}


void Foam::dgSSPRK3TimeScheme::clearResiduals()
{
    forAll(scalarStates_, i)
    {
        scalarStates_[i].clearResiduals();
    }

    forAll(vectorStates_, i)
    {
        vectorStates_[i].clearResiduals();
    }
}


void Foam::dgSSPRK3TimeScheme::updateStage
(
    const label cellID,
    const scalar dt
)
{
    validateStageIndex();

    forAll(scalarStates_, i)
    {
        updateSSPRK3Stage(scalarStates_[i], mesh_, stage(), cellID, dt);
    }

    forAll(vectorStates_, i)
    {
        updateSSPRK3Stage(vectorStates_[i], mesh_, stage(), cellID, dt);
    }
}


void Foam::dgSSPRK3TimeScheme::finalizeTimeStep
(
    const label cellID,
    const scalar dt
)
{
    // SSPRK3 completes the final update in stage 2.
}


void Foam::dgSSPRK3TimeScheme::reset()
{
    forAll(scalarStates_, i)
    {
        scalarStates_[i].reset();
    }

    forAll(vectorStates_, i)
    {
        vectorStates_[i].reset();
    }
}

} // End namespace Foam

// ************************************************************************* //
