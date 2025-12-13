/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DGFoam: Discontinuous Galerkin CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | GPU-friendly CFD solver framework
     \\/     M anipulation  |
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgHLLEFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dgHLLEFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgHLLEFluxSolver, dictionary);


// * * * * * * * * Constructors * * * * * * * * //

Foam::dgHLLEFluxSolver::dgHLLEFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgFluxSolver(name, dict, mesh),

    rho_( mesh.getFvMesh().lookupObject<dgField<scalar>>("rho") ),
    U_  ( mesh.getFvMesh().lookupObject<dgField<vector>>("U")   ),
    p_  ( mesh.getFvMesh().lookupObject<dgField<scalar>>("p")   ),
    a_  ( mesh.getFvMesh().lookupObject<dgField<scalar>>("a")   ),
    h_  ( mesh.getFvMesh().lookupObject<dgField<scalar>>("h")   ),
    gamma_( mesh.getFvMesh().lookupObject<dgField<scalar>>("gamma") )
{
    read(dict);
}


// * * * * * * * * Read dictionary * * * * * * * //

void Foam::dgHLLEFluxSolver::read(const dictionary& dict)
{
    word s = dict.lookupOrDefault<word>("speedEstimate", "davis");

    if (s == "davis")
        speedEst_ = seDavis;
    else if (s == "roeEinfeldt")
        speedEst_ = seRoeEinfeldt;
    else
    {
        WarningInFunction
            << "Unknown speedEstimate \"" << s
            << "\", using 'davis'." << nl;

        speedEst_ = seDavis;
    }
}


// * * * * * * * * Wave speed calculation * * * * * * * //

void Foam::dgHLLEFluxSolver::calcWaveSpeed
(
    scalar& SL,
    scalar& SR,

    const scalar rhoL,
    const scalar rhoR,
    const vector& ULv,
    const vector& URv,
    const scalar pL,
    const scalar pR,
    const scalar aL,
    const scalar aR,
    const scalar gammaL,
    const scalar gammaR,
    const scalar hL,
    const scalar hR,
    const vector& n
) const
{
    const scalar UnL = (ULv & n);
    const scalar UnR = (URv & n);

    if (speedEst_ == seDavis)
    {
        SL = min(UnL - aL, UnR - aR);
        SR = max(UnL + aL, UnR + aR);
        return;
    }

    // Roe-Einfeldt
    const scalar sL = sqrt(max(rhoL, SMALL));
    const scalar sR = sqrt(max(rhoR, SMALL));
    const scalar denom = sL + sR + SMALL;

    const vector Uroe = (sL*ULv + sR*URv) / denom;
    const scalar UnRoe = (Uroe & n);

    const scalar HL = hL + 0.5*magSqr(ULv);
    const scalar HR = hR + 0.5*magSqr(URv);
    const scalar Hroe = (sL*HL + sR*HR) / denom;

    const scalar gm1L = gammaL - 1.0;
    const scalar gm1R = gammaR - 1.0;

    const scalar gm1Roe =
    (
        gm1L*HL + gm1R*HR
    ) / max(HL + HR, SMALL);

    const scalar aRoe2 =
        max(gm1Roe*(Hroe - 0.5*magSqr(Uroe)), SMALL);

    const scalar aRoe = sqrt(aRoe2);

    SL = min(UnL - aL, UnRoe - aRoe);
    SR = max(UnR + aR, UnRoe + aRoe);
}


// * * * * * * * * Scalar U, vector F * * * * * * * * //

void Foam::dgHLLEFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<vector>& F,
    const faceGaussField<scalar>& U
)
{
    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
    const faceGaussField<vector>& UF   = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& pF   = p_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& aF   = a_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& hF   = h_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& gF   = gamma_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            // Left/right states
            scalar rhoL = rhoF.minusValueOnFace(fI, gI);
            scalar rhoR = rhoF.plusValueOnFace(fI, gI);

            vector ULv = UF.minusValueOnFace(fI, gI);
            vector URv = UF.plusValueOnFace(fI, gI);

            scalar pL = pF.minusValueOnFace(fI, gI);
            scalar pR = pF.plusValueOnFace(fI, gI);

            scalar aL = aF.minusValueOnFace(fI, gI);
            scalar aR = aF.plusValueOnFace(fI, gI);

            scalar hL = hF.minusValueOnFace(fI, gI);
            scalar hR = hF.plusValueOnFace(fI, gI);

            scalar gL = gF.minusValueOnFace(fI, gI);
            scalar gR = gF.plusValueOnFace(fI, gI);

            scalar SL, SR;
            calcWaveSpeed(SL, SR,
                          rhoL, rhoR, ULv, URv, pL, pR,
                          aL, aR, gL, gR, hL, hR, n);

            const vector FL = F.minusValueOnFace(fI, gI);
            const vector FR = F.plusValueOnFace(fI, gI);

            const scalar fL = (FL & n);
            const scalar fR = (FR & n);

            const scalar UL = U.minusValueOnFace(fI, gI);
            const scalar UR = U.plusValueOnFace(fI, gI);

            scalar fn;

            if (SL >= 0)
                fn = fL;
            else if (SR <= 0)
                fn = fR;
            else
            {
                fn = (SR*fL - SL*fR + SL*SR*(UR - UL))
                   / (SR - SL + VSMALL);
            }

            F.fluxOnFace(fI, gI) = fn * n;
        }
    }
}


// * * * * * * * * Vector U, tensor F * * * * * * * * //

void Foam::dgHLLEFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>& U
)
{
    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
    const faceGaussField<vector>& UF   = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& pF   = p_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& aF   = a_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& hF   = h_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& gF   = gamma_.gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            scalar rhoL = rhoF.minusValueOnFace(fI, gI);
            scalar rhoR = rhoF.plusValueOnFace(fI, gI);

            vector ULv = UF.minusValueOnFace(fI, gI);
            vector URv = UF.plusValueOnFace(fI, gI);

            scalar pL = pF.minusValueOnFace(fI, gI);
            scalar pR = pF.plusValueOnFace(fI, gI);

            scalar aL = aF.minusValueOnFace(fI, gI);
            scalar aR = aF.plusValueOnFace(fI, gI);

            scalar hL = hF.minusValueOnFace(fI, gI);
            scalar hR = hF.plusValueOnFace(fI, gI);

            scalar gL = gF.minusValueOnFace(fI, gI);
            scalar gR = gF.plusValueOnFace(fI, gI);

            scalar SL, SR;
            calcWaveSpeed(SL, SR,
                          rhoL, rhoR, ULv, URv, pL, pR,
                          aL, aR, gL, gR, hL, hR, n);

            const tensor FL = F.minusValueOnFace(fI, gI);
            const tensor FR = F.plusValueOnFace(fI, gI);

            const vector fL = (FL & n);
            const vector fR = (FR & n);

            const vector UL = U.minusValueOnFace(fI, gI);
            const vector UR = U.plusValueOnFace(fI, gI);

            vector fn;

            if (SL >= 0)
                fn = fL;
            else if (SR <= 0)
                fn = fR;
            else
            {
                fn =
                (
                    SR*fL - SL*fR + SL*SR*(UR - UL)
                ) / (SR - SL + VSMALL);
            }

            // IMPORTANT:
            // For vector U and tensor F, flux is a vector (do NOT multiply by n)
            F.fluxOnFace(fI, gI) = fn;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //

