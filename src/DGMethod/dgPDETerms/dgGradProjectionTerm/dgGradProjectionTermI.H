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

template<class UType>
inline void Foam::dgGradProjectionTerm<UType>::checkTypes() const
{
    static_assert
    (
        std::is_same<UType, scalar>::value
     || std::is_same<UType, vector>::value,
        "dgGradProjectionTerm: invalid UType. Allowed: scalar or vector"
    );

    static_assert
    (
        std::is_same
        <
            residual_type,
            typename outerProduct<vector, UType>::type
        >::value,
        "dgGradProjectionTerm: volume and surface projection result types differ"
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class UType>
Foam::dgGradProjectionTerm<UType>::dgGradProjectionTerm
(
    const word& termName,
    const dgGeomMesh& dgMesh,
    const dgFluxSolver::equationType eqnType
)
:
    termName_(termName),
    eqnType_(eqnType),
    dgMesh_(dgMesh),
    U_(nullptr)
{
    checkTypes();
}


template<class UType>
void Foam::dgGradProjectionTerm<UType>::add
(
    GaussField<UType>& U
)
{
    U_ = &U;
}

template<class UType>
void Foam::dgGradProjectionTerm<UType>::clear()
{
    U_ = nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class UType>
void Foam::dgGradProjectionTerm<UType>::calculateFlux(const label cellID)
{
    if (!U_)
    {
        FatalErrorInFunction
            << "No active state field set for " << termName_
            << ". Call add(U) before calculateFlux(cellID)."
            << abort(FatalError);
    }

    if (U_->cellID() != cellID)
    {
        FatalErrorInFunction
            << "Active state field for " << termName_ << " belongs to cell "
            << U_->cellID()
            << ", but calculateFlux was called for cell " << cellID
            << abort(FatalError);
    }

    faceGaussField<UType>& C = U_->faceField();
    C.initFlux();

    for (label faceI = 0; faceI < C.nFaces(); ++faceI)
    {
        const vector& n = C.normals()[faceI];
        const bool useBoundaryMid =
            C.isBoundary(faceI) && !C.isProcessorPatch(faceI);

        for (label g = 0; g < C.nGaussPerFace(); ++g)
        {
            const UType CStar =
                useBoundaryMid
              ? C.midValueOnFace(faceI, g)
              : scalar(0.5)
               *(
                    C.minusValueOnFace(faceI, g)
                  + C.plusValueOnFace(faceI, g)
                );

            if constexpr (std::is_same_v<UType, scalar>)
            {
                C.fluxOnFace(faceI, g) = CStar*n;
            }
            else
            {
                C.fluxOnFace(faceI, g) = CStar;
            }
        }
    }
}

// ************************************************************************* //
