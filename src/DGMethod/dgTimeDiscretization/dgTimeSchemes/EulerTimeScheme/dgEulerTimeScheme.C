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

#include "dgEulerTimeScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "pTraits.H"

namespace Foam
{

defineTypeNameAndDebug(dgEulerTimeScheme, 0);
addToRunTimeSelectionTable(dgTimeScheme, dgEulerTimeScheme, dictionary);

// * * * * * * * * * * * * * Local Helpers  * * * * * * * * * * * * * * * * //

template<class Type>
void saveOldCell(EulerFieldState<Type>& state, const label cellID)
{
    // Snapshot the current cell before assembling its update.
    dgField<Type>& U = state.U();
    dgField<Type>& UOld = state.UOld();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot save Euler state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    UOld.dof()[cellID] = U.dof()[cellID];
    UOld.gaussFields()[cellID] = U.gaussFields()[cellID];
}


template<class Type>
void restoreOldCell(EulerFieldState<Type>& state, const label cellID)
{
    // Restore one cell from the saved step-start state.
    dgField<Type>& U = state.U();
    const dgField<Type>& UOld = state.UOld();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot restore Euler state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    U.dof()[cellID] = UOld.dof()[cellID];
    U.gaussFields()[cellID] = UOld.gaussFields()[cellID];
}


template<class Type>
void clearResidualCell(EulerFieldState<Type>& state, const label cellID)
{
    // Residual is assembled fresh for each cell visit.
    List<Type>& RCell = state.R(cellID);

    forAll(RCell, dofI)
    {
        RCell[dofI] = pTraits<Type>::zero;
    }
}


template<class Type>
void updateEulerState
(
    EulerFieldState<Type>& state,
    const dgGeomMesh& mesh,
    const label cellID,
    const scalar dt
)
{
    // Explicit Euler uses only U^n and the current residual.
    dgField<Type>& U = state.U();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot update Euler state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    const dgField<Type>& UOld = state.UOld();
    const List<Type>& UOldCellDof = UOld.dof()[cellID].dof();
    const List<Type>& RCell = state.R(cellID);
    List<Type>& UCellDof = U.dof()[cellID].dof();
    const List<scalar> massDiag(mesh.cells()[cellID]->massMatrixDiag());

    forAll(UCellDof, dofI)
    {
        UCellDof[dofI] = UOldCellDof[dofI] + dt*(RCell[dofI]/massDiag[dofI]);
    }

    U.gaussFields()[cellID].interpolateFromDof();
    U.dof().updateCellDof(cellID);
}

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::dgEulerTimeScheme::dgEulerTimeScheme
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

void Foam::dgEulerTimeScheme::read(const dictionary& dict)
{
    // Euler has no extra coefficients to parse yet.
    dict_ = dict;
}

// * * * * * * * * * * * * * * * Private Helpers  * * * * * * * * * * * * * //


label Foam::dgEulerTimeScheme::findScalarIndex(const word& fieldName) const
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


label Foam::dgEulerTimeScheme::findVectorIndex(const word& fieldName) const
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


void Foam::dgEulerTimeScheme::validateStageIndex() const
{
    if (stage() != 0)
    {
        FatalErrorInFunction
            << "Euler scheme supports only stage index 0, but got stageI = "
            << stage() << "."
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Setup  * * * * * * * * * * * * * * * * * * //

void Foam::dgEulerTimeScheme::registerField(dgField<scalar>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register scalar field " << U.name()
            << " to dgEulerTimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findScalarIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Scalar field " << U.name()
            << " is already registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    const label newI = scalarNames_.size();
    scalarNames_.setSize(newI + 1);
    scalarNames_[newI] = U.name();

    scalarStates_.setSize(newI + 1);
    scalarStates_.set(newI, new EulerFieldState<scalar>(U, mesh_));
}


void Foam::dgEulerTimeScheme::registerField(dgField<vector>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register vector field " << U.name()
            << " to dgEulerTimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findVectorIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Vector field " << U.name()
            << " is already registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    const label newI = vectorNames_.size();
    vectorNames_.setSize(newI + 1);
    vectorNames_[newI] = U.name();

    vectorStates_.setSize(newI + 1);
    vectorStates_.set(newI, new EulerFieldState<vector>(U, mesh_));
}

// * * * * * * * * * * * * * * * Access  * * * * * * * * * * * * * * * * * //

label Foam::dgEulerTimeScheme::nStages() const
{
    return 1;
}


dgField<scalar>& Foam::dgEulerTimeScheme::scalarField(const word& fieldName)
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


const dgField<scalar>& Foam::dgEulerTimeScheme::scalarField(const word& fieldName) const
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


dgField<vector>& Foam::dgEulerTimeScheme::vectorField(const word& fieldName)
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


const dgField<vector>& Foam::dgEulerTimeScheme::vectorField(const word& fieldName) const
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


List<List<scalar>>& Foam::dgEulerTimeScheme::scalarResidual
(
    const word& fieldName,
    const label stageI
)
{
    validateStageIndex();

    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R();
}


const List<List<scalar>>& Foam::dgEulerTimeScheme::scalarResidual
(
    const word& fieldName,
    const label stageI
) const
{
    validateStageIndex();

    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R();
}


List<List<vector>>& Foam::dgEulerTimeScheme::vectorResidual
(
    const word& fieldName,
    const label stageI
)
{
    validateStageIndex();

    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R();
}


const List<List<vector>>& Foam::dgEulerTimeScheme::vectorResidual
(
    const word& fieldName,
    const label stageI
) const
{
    validateStageIndex();

    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgEulerTimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R();
}

// * * * * * * * * * * * * * * * Core API  * * * * * * * * * * * * * * * * //

void Foam::dgEulerTimeScheme::saveOldFields(const label cellID)
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


void Foam::dgEulerTimeScheme::restoreOldFields(const label cellID)
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


void Foam::dgEulerTimeScheme::saveStageFields(const label cellID)
{
    // Euler has no intermediate stage field.
    validateStageIndex();
}


void Foam::dgEulerTimeScheme::restoreStageFields(const label cellID)
{
    // Euler has no intermediate stage field.
    validateStageIndex();
}


void Foam::dgEulerTimeScheme::clearResiduals(const label cellID)
{
    forAll(scalarStates_, i)
    {
        clearResidualCell(scalarStates_[i], cellID);
    }

    forAll(vectorStates_, i)
    {
        clearResidualCell(vectorStates_[i], cellID);
    }
}


void Foam::dgEulerTimeScheme::clearResiduals()
{
    // Utility path for full reset over the whole mesh.
    for (label cellID = 0; cellID < mesh_.nCells(); ++cellID)
    {
        clearResiduals(cellID);
    }
}


void Foam::dgEulerTimeScheme::updateStage
(
    const label cellID,
    const scalar dt
)
{
    // One cell is advanced after its residual has been assembled.
    validateStageIndex();

    forAll(scalarStates_, i)
    {
        updateEulerState(scalarStates_[i], mesh_, cellID, dt);
    }

    forAll(vectorStates_, i)
    {
        updateEulerState(vectorStates_[i], mesh_, cellID, dt);
    }
}


void Foam::dgEulerTimeScheme::finalizeTimeStep
(
    const label cellID,
    const scalar dt
)
{
    // Euler completes the time-step update in its only stage.
}


void Foam::dgEulerTimeScheme::reset()
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
