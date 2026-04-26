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
    const dgBasisField& basis,
    const GaussField<UType>& U,
    GaussField<FType>& F,
    dgFluxSolverManager& fluxSolver,
    const bool removeHighOrder
)
:
    termName_(termName),
    eqnType_(eqnType),
    fluxType_(fluxType),
    cellID_(basis.cellID()),
    basis_(basis),
    U_(U),
    F_(new GaussField<FType>(F)),
    fluxSolver_(fluxSolver),
    nDof_(basis_.nDof()),
    removeHighOrder_(removeHighOrder),
    R_(nDof_, pTraits<UType>::zero)
{
    checkTypes();

    F_.ref().faceField().initFlux();

    fluxSolver_.solver(termName_+"Flux", fluxType_).setEquationType(eqnType_);
}


template<class UType, class FType>
Foam::dgDiffusiveTerm<UType, FType>::dgDiffusiveTerm
(
    const word& termName,
    dgFluxSolver::equationType eqnType,
    dgFluxSolver::fluxType fluxType,
    const dgBasisField& basis,
    const GaussField<UType>& U,
    const tmp<GaussField<FType>>& tF,
    dgFluxSolverManager& fluxSolver,
    const bool removeHighOrder
)
:
    termName_(termName),
    eqnType_(eqnType),
    fluxType_(fluxType),
    cellID_(basis.cellID()),
    basis_(basis),
    U_(U),
    F_(tF),
    fluxSolver_(fluxSolver),
    nDof_(basis_.nDof()),
    removeHighOrder_(removeHighOrder),
    R_(nDof_, pTraits<UType>::zero)
{
    checkTypes();

    F_.ref().faceField().initFlux();

    fluxSolver_.solver(termName_+"Flux", fluxType_).setEquationType(eqnType_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::updateFlux()
{
    fluxSolver_.solver(termName_+"Flux", fluxType_).computeFlux
    (
        cellID_,
        F_.ref().faceField(),
        U_.faceField()
    );
}


template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::calcResidual()
{
    forAll (R_, i)
    {
        R_[i] = pTraits<UType>::zero;
    }

    const label nResidualDof = removeHighOrder_ ? min(nDof_, label(1)) : nDof_;

    for (label i = 0; i < nResidualDof; ++i)
    {
        UType volTerm = pTraits<UType>::zero;

        if (i != 0)
        {
            const GaussField<vector>& gradPhi = basis_.getDBasis(i);
            volTerm = Foam::volumeIntegral
            (
                gradPhi.cellField(),
                F_().cellField()
            );
        }

        const GaussField<scalar>& phi = basis_.getBasis(i);

        const UType surfTerm = Foam::surfaceIntegral(phi.faceField(), F_().faceField());

        R_[i] += volTerm;
        R_[i] -= surfTerm;
    }
}

// ************************************************************************* //
