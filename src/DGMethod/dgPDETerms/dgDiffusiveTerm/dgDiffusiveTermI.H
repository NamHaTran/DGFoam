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


template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::updateCentralFlux
(
    faceGaussField<vector>& F
)
{
    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            if (F.isBoundary(fI) && !F.isProcessorPatch(fI))
            {
                F.plusValueOnFace(fI, gI) = F.minusValueOnFace(fI, gI);
            }

            const vector FL = F.minusValueOnFace(fI, gI);
            const vector FR = F.plusValueOnFace(fI, gI);
            const scalar fC = scalar(0.5)*((FL & n) + (FR & n));

            F.fluxOnFace(fI, gI) = fC*n;
        }
    }
}


template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::updateCentralFlux
(
    faceGaussField<tensor>& F
)
{
    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            if (F.isBoundary(fI) && !F.isProcessorPatch(fI))
            {
                F.plusValueOnFace(fI, gI) = F.minusValueOnFace(fI, gI);
            }

            const tensor FL = F.minusValueOnFace(fI, gI);
            const tensor FR = F.plusValueOnFace(fI, gI);
            const vector fC = scalar(0.5)*((FL & n) + (FR & n));

            F.fluxOnFace(fI, gI) = fC;
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class UType, class FType>
Foam::dgDiffusiveTerm<UType, FType>::dgDiffusiveTerm
(
    const word& termName,
    dgFluxSolver::equationType eqnType,
    const dgBasisField& basis,
    const GaussField<UType>& U,
    GaussField<FType>& F
)
:
    termName_(termName),
    eqnType_(eqnType),
    cellID_(basis.cellID()),
    basis_(basis),
    U_(U),
    F_(new GaussField<FType>(F)),
    nDof_(basis_.nDof()),
    R_(nDof_, pTraits<UType>::zero)
{
    checkTypes();

    F_.ref().faceField().initFlux();
}


template<class UType, class FType>
Foam::dgDiffusiveTerm<UType, FType>::dgDiffusiveTerm
(
    const word& termName,
    dgFluxSolver::equationType eqnType,
    const dgBasisField& basis,
    const GaussField<UType>& U,
    const tmp<GaussField<FType>>& tF
)
:
    termName_(termName),
    eqnType_(eqnType),
    cellID_(basis.cellID()),
    basis_(basis),
    U_(U),
    F_(tF),
    nDof_(basis_.nDof()),
    R_(nDof_, pTraits<UType>::zero)
{
    checkTypes();

    F_.ref().faceField().initFlux();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::updateFlux()
{
    updateCentralFlux(F_.ref().faceField());
}


template<class UType, class FType>
void Foam::dgDiffusiveTerm<UType, FType>::calcResidual()
{
    forAll (R_, i)
    {
        R_[i] = pTraits<UType>::zero;
    }

    for (label i = 0; i < nDof_; ++i)
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
