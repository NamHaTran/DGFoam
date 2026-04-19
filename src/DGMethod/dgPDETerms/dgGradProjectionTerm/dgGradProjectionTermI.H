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
    const dgBasisField& basis,
    const GaussField<UType>& U,
    const dgFluxSolver::equationType eqnType
)
:
    termName_(termName),
    eqnType_(eqnType),
    cellID_(basis.cellID()),
    basis_(basis),
    U_(new GaussField<UType>(U)),
    nDof_(basis_.nDof()),
    R_(nDof_, pTraits<residual_type>::zero)
{
    checkTypes();
}


template<class UType>
Foam::dgGradProjectionTerm<UType>::dgGradProjectionTerm
(
    const word& termName,
    const dgBasisField& basis,
    const tmp<GaussField<UType>>& tU,
    const dgFluxSolver::equationType eqnType
)
:
    termName_(termName),
    eqnType_(eqnType),
    cellID_(basis.cellID()),
    basis_(basis),
    U_(tU),
    nDof_(basis_.nDof()),
    R_(nDof_, pTraits<residual_type>::zero)
{
    checkTypes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class UType>
void Foam::dgGradProjectionTerm<UType>::updateFlux()
{
    faceGaussField<UType>& UFace = U_.ref().faceField();

    const label nFaces = UFace.nFaces();
    const label nGauss = UFace.nGaussPerFace();

    for (label faceI = 0; faceI < nFaces; ++faceI)
    {
        if (UFace.isBoundary(faceI) && !UFace.isProcessorPatch(faceI))
        {
            continue;
        }

        for (label gI = 0; gI < nGauss; ++gI)
        {
            UFace.midValueOnFace(faceI, gI) =
                scalar(0.5)
               *(
                    UFace.minusValueOnFace(faceI, gI)
                  + UFace.plusValueOnFace(faceI, gI)
                );
        }
    }
}


template<class UType>
void Foam::dgGradProjectionTerm<UType>::calcResidual()
{
    forAll (R_, i)
    {
        R_[i] = pTraits<residual_type>::zero;
    }

    for (label i = 0; i < nDof_; ++i)
    {
        residual_type volTerm = pTraits<residual_type>::zero;

        if (i != 0)
        {
            const GaussField<vector>& gradPhi = basis_.getDBasis(i);
            volTerm = Foam::gradProjectionVolumeIntegral
            (
                gradPhi.cellField(),
                U_().cellField()
            );
        }

        const GaussField<scalar>& phi = basis_.getBasis(i);

        const residual_type surfTerm =
            Foam::gradProjectionSurfaceIntegral
            (
                phi.faceField(),
                U_().faceField()
            );

        R_[i] -= volTerm;
        R_[i] += surfTerm;
    }
}

// ************************************************************************* //
