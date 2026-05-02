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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class UType, class FType>
inline void Foam::dgDiffusiveTerm<UType, FType>::checkTypes() const
{
    static_assert
    (
        (std::is_same<UType, scalar>::value && std::is_same<FType, vector>::value)
     ||
        (std::is_same<UType, vector>::value && std::is_same<FType, tensor>::value),
        "dgDiffusiveTerm: invalid (UType, FType) combination. "
        "Allowed: (scalar, vector) or (vector, tensor)"
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class UType, class FType>
Foam::dgDiffusiveTerm<UType, FType>::dgDiffusiveTerm
(
    const word& termName,
    dgFluxSolver::equationType eqnType,
    dgFluxSolver::fluxType fluxType,
    const dgGeomMesh& dgMesh,
    dgFluxSolverManager& fluxSolver
)
:
    termName_(termName),
    eqnType_(eqnType),
    fluxType_(fluxType),
    dgMesh_(dgMesh),
    U_(nullptr),
    F_(nullptr),
    fluxSolver_(fluxSolver),
    termFluxSolver_(fluxSolver_.solver(termName_+"Flux", fluxType_))
{
    checkTypes();

    termFluxSolver_.setEquationType(eqnType_);
}


template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::add
(
    const GaussField<UType>& U,
    GaussField<FType>& F
)
{
    U_ = &U;
    F_ = &F;
    F_->faceField().initFlux();
}

template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::add
(
    const GaussField<UType>& U,
    tmp<GaussField<FType>>& tF
)
{
    add(U, tF.ref());
}

template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::clear()
{
    U_ = nullptr;
    F_ = nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::updateFlux(const label cellID)
{
    if (!U_ || !F_)
    {
        FatalErrorInFunction
            << "No active fields set for " << termName_
            << ". Call add(U, F) before updateFlux(cellID)."
            << abort(FatalError);
    }

    if (U_->cellID() != cellID || F_->cellID() != cellID)
    {
        FatalErrorInFunction
            << "Active fields for " << termName_ << " belong to cell "
            << U_->cellID() << "/" << F_->cellID()
            << ", but updateFlux was called for cell " << cellID
            << abort(FatalError);
    }

    termFluxSolver_.computeFlux
    (
        cellID,
        F_->faceField(),
        U_->faceField()
    );
}


// ************************************************************************* //
