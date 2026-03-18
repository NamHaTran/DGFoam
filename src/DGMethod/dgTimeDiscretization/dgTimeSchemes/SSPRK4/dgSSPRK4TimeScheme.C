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

#include "dgSSPRK4TimeScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "pTraits.H"

namespace Foam
{

defineTypeNameAndDebug(dgSSPRK4TimeScheme, 0);
addToRunTimeSelectionTable(dgTimeScheme, dgSSPRK4TimeScheme, dictionary);

// * * * * * * * * * * * * * Local Helpers  * * * * * * * * * * * * * * * * //

template<class Type>
void copyCell(dgField<Type>& dst, const dgField<Type>& src, const label cellID)
{
    // Cell-local copy keeps the scheme inside the solver cell loop.
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
void syncUpdatedCell(dgField<Type>& U, const label cellID)
{
    // Update the Gauss representation after changing the DoFs.
    U.gaussFields()[cellID].interpolateFromDof();
    U.dof().updateCellDof(cellID);
}


template<class Type>
void saveOldCell(SSPRK4FieldState<Type>& state, const label cellID)
{
    copyCell(state.UOld(), state.U(), cellID);
}


template<class Type>
void restoreOldCell(SSPRK4FieldState<Type>& state, const label cellID)
{
    copyCell(state.U(), state.UOld(), cellID);
}


template<class Type>
void saveStageCell
(
    SSPRK4FieldState<Type>& state,
    const label stageFieldI,
    const label cellID
)
{
    copyCell(state.UStage(stageFieldI), state.U(), cellID);
}


template<class Type>
void restoreStageCell
(
    SSPRK4FieldState<Type>& state,
    const label stageFieldI,
    const label cellID
)
{
    copyCell(state.U(), state.UStage(stageFieldI), cellID);
}


template<class Type>
void clearResidualCell
(
    SSPRK4FieldState<Type>& state,
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
void updateSSPRK4Stage
(
    SSPRK4FieldState<Type>& state,
    const dgGeomMesh& mesh,
    const label stageI,
    const label cellID,
    const scalar dt
)
{
    // Each stage uses the stored SSPRK4 intermediates for one cell only.
    dgField<Type>& U = state.U();

    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot update SSPRK4 state for field " << U.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    List<Type>& UCellDof = U.dof()[cellID].dof();
    const List<scalar> massDiag(mesh.cells()[cellID]->massMatrixDiag());

    switch (stageI)
    {
        case 0:
        {
            const List<Type>& U0 = state.UOld().dof()[cellID].dof();
            const List<Type>& R0 = state.R(0, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    U0[dofI] + 0.391752226571890*dt*(R0[dofI]/massDiag[dofI]);
            }

            syncUpdatedCell(U, cellID);
            copyCell(state.UStage(0), U, cellID);
            break;
        }

        case 1:
        {
            const List<Type>& U0 = state.UOld().dof()[cellID].dof();
            const List<Type>& U1 = state.UStage(0).dof()[cellID].dof();
            const List<Type>& R1 = state.R(1, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    0.444370493651235*U0[dofI]
                  + 0.555629506348765*U1[dofI]
                  + 0.368410593050371*dt*(R1[dofI]/massDiag[dofI]);
            }

            syncUpdatedCell(U, cellID);
            copyCell(state.UStage(1), U, cellID);
            break;
        }

        case 2:
        {
            const List<Type>& U0 = state.UOld().dof()[cellID].dof();
            const List<Type>& U2 = state.UStage(1).dof()[cellID].dof();
            const List<Type>& R2 = state.R(2, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    0.620101851488403*U0[dofI]
                  + 0.379898148511597*U2[dofI]
                  + 0.251891774247377*dt*(R2[dofI]/massDiag[dofI]);
            }

            syncUpdatedCell(U, cellID);
            copyCell(state.UStage(2), U, cellID);
            break;
        }

        case 3:
        {
            const List<Type>& U0 = state.UOld().dof()[cellID].dof();
            const List<Type>& U3 = state.UStage(2).dof()[cellID].dof();
            const List<Type>& R3 = state.R(3, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    0.178079954393132*U0[dofI]
                  + 0.821920045606868*U3[dofI]
                  + 0.544974750228521*dt*(R3[dofI]/massDiag[dofI]);
            }

            syncUpdatedCell(U, cellID);
            copyCell(state.UStage(3), U, cellID);
            break;
        }

        case 4:
        {
            // Final stage combines earlier stored states and the last residual.
            const List<Type>& U2 = state.UStage(1).dof()[cellID].dof();
            const List<Type>& U3 = state.UStage(2).dof()[cellID].dof();
            const List<Type>& U4 = state.UStage(3).dof()[cellID].dof();
            const List<Type>& R3 = state.R(3, cellID);
            const List<Type>& R4 = state.R(4, cellID);

            forAll(UCellDof, dofI)
            {
                UCellDof[dofI] =
                    0.517231671970585*U2[dofI]
                  + 0.096059710526147*U3[dofI]
                  + 0.063692468666290*dt*(R3[dofI]/massDiag[dofI])
                  + 0.386708617503269*U4[dofI]
                  + 0.226007483236906*dt*(R4[dofI]/massDiag[dofI]);
            }

            syncUpdatedCell(U, cellID);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Invalid SSPRK4 stage index " << stageI
                << ". Valid stage indices are 0..4."
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::dgSSPRK4TimeScheme::dgSSPRK4TimeScheme
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

void Foam::dgSSPRK4TimeScheme::read(const dictionary& dict)
{
    // SSPRK4 currently has no extra dictionary coefficients.
    dict_ = dict;
}


label Foam::dgSSPRK4TimeScheme::findScalarIndex(const word& fieldName) const
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


label Foam::dgSSPRK4TimeScheme::findVectorIndex(const word& fieldName) const
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


void Foam::dgSSPRK4TimeScheme::validateStageIndex(const label stageI) const
{
    if (stageI < 0 || stageI >= nStages())
    {
        FatalErrorInFunction
            << "SSPRK4 scheme supports only stage indices 0.." << (nStages() - 1)
            << ", but got stageI = " << stageI << "."
            << abort(FatalError);
    }
}


void Foam::dgSSPRK4TimeScheme::validateStageIndex() const
{
    validateStageIndex(stage());
}


void Foam::dgSSPRK4TimeScheme::registerField(dgField<scalar>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register scalar field " << U.name()
            << " to dgSSPRK4TimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findScalarIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Scalar field " << U.name()
            << " is already registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    const label newI = scalarNames_.size();
    scalarNames_.setSize(newI + 1);
    scalarNames_[newI] = U.name();
    scalarStates_.setSize(newI + 1);
    scalarStates_.set(newI, new SSPRK4FieldState<scalar>(U, mesh_));
}


void Foam::dgSSPRK4TimeScheme::registerField(dgField<vector>& U)
{
    if (!U.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register vector field " << U.name()
            << " to dgSSPRK4TimeScheme because it has no DOF storage."
            << abort(FatalError);
    }

    if (findVectorIndex(U.name()) >= 0)
    {
        FatalErrorInFunction
            << "Vector field " << U.name()
            << " is already registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    const label newI = vectorNames_.size();
    vectorNames_.setSize(newI + 1);
    vectorNames_[newI] = U.name();
    vectorStates_.setSize(newI + 1);
    vectorStates_.set(newI, new SSPRK4FieldState<vector>(U, mesh_));
}


label Foam::dgSSPRK4TimeScheme::nStages() const
{
    return 5;
}


dgField<scalar>& Foam::dgSSPRK4TimeScheme::scalarField(const word& fieldName)
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


const dgField<scalar>& Foam::dgSSPRK4TimeScheme::scalarField(const word& fieldName) const
{
    const label i = findScalarIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].U();
}


dgField<vector>& Foam::dgSSPRK4TimeScheme::vectorField(const word& fieldName)
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


const dgField<vector>& Foam::dgSSPRK4TimeScheme::vectorField(const word& fieldName) const
{
    const label i = findVectorIndex(fieldName);

    if (i < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].U();
}


List<List<scalar>>& Foam::dgSSPRK4TimeScheme::scalarResidual
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
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R(stageI);
}


const List<List<scalar>>& Foam::dgSSPRK4TimeScheme::scalarResidual
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
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return scalarStates_[i].R(stageI);
}


List<List<vector>>& Foam::dgSSPRK4TimeScheme::vectorResidual
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
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R(stageI);
}


const List<List<vector>>& Foam::dgSSPRK4TimeScheme::vectorResidual
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
            << " is not registered in dgSSPRK4TimeScheme."
            << abort(FatalError);
    }

    return vectorStates_[i].R(stageI);
}


void Foam::dgSSPRK4TimeScheme::saveOldFields(const label cellID)
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


void Foam::dgSSPRK4TimeScheme::restoreOldFields(const label cellID)
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


void Foam::dgSSPRK4TimeScheme::saveStageFields(const label cellID)
{
    // Stages 0..3 are stored because later stages reuse them.
    validateStageIndex();

    if (stage() < 4)
    {
        forAll(scalarStates_, i)
        {
            saveStageCell(scalarStates_[i], stage(), cellID);
        }

        forAll(vectorStates_, i)
        {
            saveStageCell(vectorStates_[i], stage(), cellID);
        }
    }
}


void Foam::dgSSPRK4TimeScheme::restoreStageFields(const label cellID)
{
    // This is mainly useful when replaying a stored intermediate stage.
    validateStageIndex();

    if (stage() < 4)
    {
        forAll(scalarStates_, i)
        {
            restoreStageCell(scalarStates_[i], stage(), cellID);
        }

        forAll(vectorStates_, i)
        {
            restoreStageCell(vectorStates_[i], stage(), cellID);
        }
    }
}


void Foam::dgSSPRK4TimeScheme::clearResiduals(const label cellID)
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


void Foam::dgSSPRK4TimeScheme::clearResiduals()
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


void Foam::dgSSPRK4TimeScheme::updateStage
(
    const label cellID,
    const scalar dt
)
{
    validateStageIndex();

    forAll(scalarStates_, i)
    {
        updateSSPRK4Stage(scalarStates_[i], mesh_, stage(), cellID, dt);
    }

    forAll(vectorStates_, i)
    {
        updateSSPRK4Stage(vectorStates_[i], mesh_, stage(), cellID, dt);
    }
}


void Foam::dgSSPRK4TimeScheme::finalizeTimeStep
(
    const label cellID,
    const scalar dt
)
{
    // SSPRK4 completes the final update in stage 4.
}


void Foam::dgSSPRK4TimeScheme::reset()
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
