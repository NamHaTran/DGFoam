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

#include "dgRK4TimeScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "pTraits.H"

namespace Foam
{

defineTypeNameAndDebug(dgRK4TimeScheme, 0);
addToRunTimeSelectionTable(dgTimeScheme, dgRK4TimeScheme, dictionary);

// * * * * * * * * * * * * * Local Helpers  * * * * * * * * * * * * * * * * //

template<class Type>
void saveOldCell(RK4FieldState<Type>& state, const label cellID)
{
    // Keep a per-cell copy of U^n for all later RK combinations.
    dgField<Type>& U = state.U();
    dgField<Type>& UOld = state.UOld();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot save RK4 state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    UOld.dof()[cellID] = U.dof()[cellID];
    UOld.gaussFields()[cellID] = U.gaussFields()[cellID];
}


template<class Type>
void restoreOldCell(RK4FieldState<Type>& state, const label cellID)
{
    // Restore one cell from the saved step-start state.
    dgField<Type>& U = state.U();
    const dgField<Type>& UOld = state.UOld();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot restore RK4 state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    U.dof()[cellID] = UOld.dof()[cellID];
    U.gaussFields()[cellID] = UOld.gaussFields()[cellID];
}


template<class Type>
void clearResidualCell
(
    RK4FieldState<Type>& state,
    const label stageI,
    const label cellID
)
{
    // Only the active stage residual is cleared here.
    List<Type>& RCell = state.R(stageI, cellID);

    forAll(RCell, dofI)
    {
        RCell[dofI] = pTraits<Type>::zero;
    }
}


template<class Type>
void updateRK4Stage
(
    RK4FieldState<Type>& state,
    const dgGeomMesh& mesh,
    const label stageI,
    const label cellID,
    const scalar dt
)
{
    // Rebuild the stage state directly from U^n and stored residuals.
    dgField<Type>& U = state.U();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot update RK4 state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    const dgField<Type>& UOld = state.UOld();
    const List<Type>& UOldCellDof = UOld.dof()[cellID].dof();
    List<Type>& UCellDof = U.dof()[cellID].dof();
    const List<scalar>& massDiag = mesh.cells()[cellID]->massMatrixDiag();

    switch (stageI)
    {
        case 0:
        {
            const List<Type>& R1 = state.R(0, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    UOldCellDof[dofI] + 0.5*dt*(R1[dofI]/massDiag[dofI]);
            }
            break;
        }

        case 1:
        {
            const List<Type>& R2 = state.R(1, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    UOldCellDof[dofI] + 0.5*dt*(R2[dofI]/massDiag[dofI]);
            }
            break;
        }

        case 2:
        {
            const List<Type>& R3 = state.R(2, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    UOldCellDof[dofI] + dt*(R3[dofI]/massDiag[dofI]);
            }
            break;
        }

        case 3:
        {
            // The last RK4 stage also performs the final weighted combine.
            const List<Type>& R1 = state.R(0, cellID);
            const List<Type>& R2 = state.R(1, cellID);
            const List<Type>& R3 = state.R(2, cellID);
            const List<Type>& R4 = state.R(3, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    UOldCellDof[dofI]
                  + (dt/6.0)
                   *
                    (
                        (R1[dofI]/massDiag[dofI])
                      + 2.0*(R2[dofI]/massDiag[dofI])
                      + 2.0*(R3[dofI]/massDiag[dofI])
                      + (R4[dofI]/massDiag[dofI])
                    );
            }
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Invalid RK4 stage index " << stageI
                << ". Valid stage indices are 0..3."
                << abort(FatalError);
        }
    }

    U.dof().updateCellDof(cellID);
}

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::dgRK4TimeScheme::dgRK4TimeScheme
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dgRK4TimeScheme::read(const dictionary& dict)
{
    // RK4 has no extra dictionary entries yet.
    dict_ = dict;
}

// * * * * * * * * * * * * * * * Private Helpers  * * * * * * * * * * * * * //

label Foam::dgRK4TimeScheme::findScalarIndex(const word& fieldName) const
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


label Foam::dgRK4TimeScheme::findVectorIndex(const word& fieldName) const
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


void Foam::dgRK4TimeScheme::validateStageIndex(const label stageI) const
{
    if (stageI < 0 || stageI >= nStages())
    {
        FatalErrorInFunction
            << "RK4 scheme supports only stage indices 0.." << (nStages() - 1)
            << ", but got stageI = " << stageI << "."
            << abort(FatalError);
    }
}


void Foam::dgRK4TimeScheme::validateStageIndex() const
{
    validateStageIndex(stage());
}

// * * * * * * * * * * * * * * * Setup  * * * * * * * * * * * * * * * * * * //

void Foam::dgRK4TimeScheme::registerField(dgField<scalar>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register scalar field " << U.name()
            << " to dgRK4TimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findScalarIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Scalar field " << U.name()
            << " is already registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    const label newI = scalarNames_.size();
    scalarNames_.setSize(newI + 1);
    scalarNames_[newI] = U.name();

    scalarStates_.setSize(newI + 1);
    scalarStates_.set(newI, new RK4FieldState<scalar>(U, mesh_));
}


void Foam::dgRK4TimeScheme::registerField(dgField<vector>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register vector field " << U.name()
            << " to dgRK4TimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findVectorIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Vector field " << U.name()
            << " is already registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    const label newI = vectorNames_.size();
    vectorNames_.setSize(newI + 1);
    vectorNames_[newI] = U.name();

    vectorStates_.setSize(newI + 1);
    vectorStates_.set(newI, new RK4FieldState<vector>(U, mesh_));
}

// * * * * * * * * * * * * * * * Access  * * * * * * * * * * * * * * * * * //

label Foam::dgRK4TimeScheme::nStages() const
{
    return 4;
}


dgField<scalar>& Foam::dgRK4TimeScheme::scalarField(const word& fieldName)
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


const dgField<scalar>& Foam::dgRK4TimeScheme::scalarField(const word& fieldName) const
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


dgField<vector>& Foam::dgRK4TimeScheme::vectorField(const word& fieldName)
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


const dgField<vector>& Foam::dgRK4TimeScheme::vectorField(const word& fieldName) const
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


List<List<scalar>>& Foam::dgRK4TimeScheme::scalarResidual
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
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R(stageI);
}


const List<List<scalar>>& Foam::dgRK4TimeScheme::scalarResidual
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
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R(stageI);
}


List<List<vector>>& Foam::dgRK4TimeScheme::vectorResidual
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
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R(stageI);
}


const List<List<vector>>& Foam::dgRK4TimeScheme::vectorResidual
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
            << " is not registered in dgRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R(stageI);
}

// * * * * * * * * * * * * * * * Core API  * * * * * * * * * * * * * * * * //

void Foam::dgRK4TimeScheme::saveOldFields(const label cellID)
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


void Foam::dgRK4TimeScheme::restoreOldFields(const label cellID)
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


void Foam::dgRK4TimeScheme::saveStageFields(const label cellID)
{
    // Classical RK4 rebuilds stages from UOld and residuals only.
    validateStageIndex();
}


void Foam::dgRK4TimeScheme::restoreStageFields(const label cellID)
{
    // Classical RK4 rebuilds stages from UOld and residuals only.
    validateStageIndex();
}


void Foam::dgRK4TimeScheme::clearResiduals(const label cellID)
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


void Foam::dgRK4TimeScheme::clearResiduals()
{
    // Utility path for a full scheme reset.
    forAll(scalarStates_, i)
    {
        scalarStates_[i].clearResiduals();
    }

    forAll(vectorStates_, i)
    {
        vectorStates_[i].clearResiduals();
    }
}


void Foam::dgRK4TimeScheme::updateStage
(
    const label cellID,
    const scalar dt
)
{
    validateStageIndex();

    forAll(scalarStates_, i)
    {
        updateRK4Stage(scalarStates_[i], mesh_, stage(), cellID, dt);
    }

    forAll(vectorStates_, i)
    {
        updateRK4Stage(vectorStates_[i], mesh_, stage(), cellID, dt);
    }
}


void Foam::dgRK4TimeScheme::finalizeTimeStep
(
    const label cellID,
    const scalar dt
)
{
    // Classical RK4 completes the final combination in stage 3.
}


void Foam::dgRK4TimeScheme::reset()
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
