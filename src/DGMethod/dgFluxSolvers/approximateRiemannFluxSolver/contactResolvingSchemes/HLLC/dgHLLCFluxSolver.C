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

#include "dgHLLCFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dgHLLCFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgHLLCFluxSolver, dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
Foam::dgHLLCFluxSolver::dgHLLCFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgFluxSolver(name, dict, mesh),
    rho_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("rho")
    ),
    U_
    (
        mesh.getFvMesh().lookupObject<dgField<vector>>("U")
    ),
    p_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("p")
    ),
    a_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("a")
    ),
    h_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("h")
    ),
    gamma_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("gamma")
    )
{
    read(dict);

    // Resize intermediate lists
    const label nFaces = mesh_.nFaces();

    SL_list_.resize(nFaces);
    SR_list_.resize(nFaces);
    SStar_list_.resize(nFaces);
    CL_list_.resize(nFaces);
    CR_list_.resize(nFaces);
    isStateComputed_.resize(nFaces);

    // Access geometric faces from dgGeomMesh
    const List<dgGeomFace*>& gFaces = mesh.faces();

    // Loop over each face
    for (label fI = 0; fI < nFaces; ++fI)
    {
        isStateComputed_[fI] = false;

        // Pointer to geometric face
        const dgGeomFace* gf = gFaces[fI];

        // Number of Gauss points on this face (owner-side Gauss points)
        const label nGauss = gf->gaussPointsOwner().size();

        // Resize the second dimension for this face
        SL_list_[fI].resize(nGauss);
        SR_list_[fI].resize(nGauss);
        SStar_list_[fI].resize(nGauss);
        CL_list_[fI].resize(nGauss);
        CR_list_[fI].resize(nGauss);

        // Initialize
        for (label gI = 0; gI < nGauss; ++gI)
        {
            SL_list_[fI][gI]    = 0.0;
            SR_list_[fI][gI]    = 0.0;
            SStar_list_[fI][gI] = 0.0;
            CL_list_[fI][gI]    = 0.0;
            CR_list_[fI][gI]    = 0.0;
        }
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dgHLLCFluxSolver::read(const dictionary& dict)
{
    word s = dict.lookupOrDefault<word>("speedEstimate", "davis");

    if (s == "davis")
    {
        speedEst_ = seDavis;
    }
    else if (s == "roeEinfeldt")
    {
        speedEst_ = seRoeEinfeldt;
    }
    else
    {
        WarningInFunction << "Unknown speedEstimate \"" << s
                      << "\". Fallback to 'davis'." << nl;
        speedEst_ = seDavis;
    }   
}

// * * * * * * * * * * * * * * scalar U, vector F  * * * * * * * * * * * * * //

void Foam::dgHLLCFluxSolver::calcIntermediateState
(
    const label cellID,
    const label localFaceID,
    const label localGaussID,
    const vector& n,

    // Output
    scalar& SL,
    scalar& SR,
    scalar& SStar,
    scalar& CL,
    scalar& CR
)
{
    // Access left and right side values at this Gauss point
    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();

    // Get global face ID from any faceGaussField
    const label globalFaceID = rhoF.globalFaceID(localFaceID);

    // Check whether the cell is owner, and state hasn't been calculated
    if (rhoF.isOwner(localFaceID, cellID) && !isStateComputed_[globalFaceID])
    {
        // Set flag
        isStateComputed_[globalFaceID] = true;

        // Calculate and cache intermediate states if the cell is owner
        const faceGaussField<vector>& UF   = U_.gaussFields()[cellID].faceField();
        const faceGaussField<scalar>& pF   = p_.gaussFields()[cellID].faceField();
        const faceGaussField<scalar>& aF   = a_.gaussFields()[cellID].faceField();
        const faceGaussField<scalar>& hF   = h_.gaussFields()[cellID].faceField();
        const faceGaussField<scalar>& gammaF = gamma_.gaussFields()[cellID].faceField();

        // Left (-) and right (+) states
        const scalar rhoR = rhoF.plusValueOnFace(localFaceID, localGaussID);
        const scalar rhoL = rhoF.minusValueOnFace(localFaceID, localGaussID);

        const vector URv = UF.plusValueOnFace(localFaceID, localGaussID);
        const vector ULv = UF.minusValueOnFace(localFaceID, localGaussID);

        const scalar pR = pF.plusValueOnFace(localFaceID, localGaussID);
        const scalar pL = pF.minusValueOnFace(localFaceID, localGaussID);

        const scalar aR = aF.plusValueOnFace(localFaceID, localGaussID);
        const scalar aL = aF.minusValueOnFace(localFaceID, localGaussID);

        const scalar hR = hF.plusValueOnFace(localFaceID, localGaussID);
        const scalar hL = hF.minusValueOnFace(localFaceID, localGaussID);

        const scalar gammaR = gammaF.plusValueOnFace(localFaceID, localGaussID);
        const scalar gammaL = gammaF.minusValueOnFace(localFaceID, localGaussID);

        // Normal velocities
        const scalar UnR = (URv & n);
        const scalar UnL = (ULv & n);

        if (speedEst_ == seDavis)
        {
            // Davis/Einfeldt (robust, positivity-friendly)
            SL = min(UnL - aL, UnR - aR);
            SR = max(UnL + aL, UnR + aR);
        }
        else // seRoeEinfeldt
        {
            // Roe averages for U, H, a (less diffusive)
            const scalar sL = sqrt(max(rhoL, SMALL));
            const scalar sR = sqrt(max(rhoR, SMALL));
            const scalar denom = sL + sR + SMALL;

            // Roe-averaged velocity
            const vector Uroe = (sL*ULv + sR*URv)/denom;
            const scalar UnRoe = (Uroe & n);

            // Left/right total enthalpy H = h + 0.5|U|^2
            const scalar H_L = hL + 0.5*magSqr(ULv);
            const scalar H_R = hR + 0.5*magSqr(URv);

            // Roe-averaged enthalpy
            const scalar H_Roe = (sL*H_L + sR*H_R)/denom;

            // Roe-averaged (gamma - 1) using enthalpy weights
            const scalar gm1L = gammaL - 1.0;
            const scalar gm1R = gammaR - 1.0;

            const scalar gm1Roe =
            (
                gm1L*H_L + gm1R*H_R
            ) / max(H_L + H_R, SMALL);

            // Roe sound speed: a^2 = (gammaRoe - 1)*(H_Roe - 0.5|Uroe|^2)
            const scalar Uroe2 = magSqr(Uroe);

            const scalar aRoe2 =
                max(gm1Roe*(H_Roe - 0.5*Uroe2), SMALL);

            const scalar aRoe = sqrt(aRoe2);

            // Roe-Einfeldt wave speeds
            SL = min(UnL - aL, UnRoe - aRoe);
            SR = max(UnR + aR, UnRoe + aRoe);
        }

        // Compute S* (contact wave speed)
        SStar =
        (
            pR - pL
        + rhoL*UnL*(SL - UnL)
        - rhoR*UnR*(SR - UnR)
        ) / (rhoL*(SL - UnL) - rhoR*(SR - UnR) + VSMALL);

        // Compute Ck = rhoK*(Sk - uK)/(Sk - S*)
        CL = rhoL*(SL - UnL)/(SL - SStar + VSMALL);
        CR = rhoR*(SR - UnR)/(SR - SStar + VSMALL);

        // Cache states
        SL_list_[globalFaceID][localGaussID] = SL;
        SR_list_[globalFaceID][localGaussID] = SR;
        SStar_list_[globalFaceID][localGaussID] = SStar;
        CL_list_[globalFaceID][localGaussID] = CL;
        CR_list_[globalFaceID][localGaussID] = CR;
    }
    else
    {
        // Get cached data
        // Sign is confirmmed
        SR = -SL_list_[globalFaceID][localGaussID];
        SL = -SR_list_[globalFaceID][localGaussID];
        SStar = -SStar_list_[globalFaceID][localGaussID];
        CR = CL_list_[globalFaceID][localGaussID];
        CL = CR_list_[globalFaceID][localGaussID];
    }
}

void Foam::dgHLLCFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<vector>& F,
    const faceGaussField<scalar>& U
)
{
    scalar SL(0.0);
    scalar SR(0.0);
    scalar SStar(0.0);
    scalar CL(0.0);
    scalar CR(0.0);

    const label nFaces(F.nFaces());
    const label nGaussPerFace(F.nGaussPerFace());

    // Loop over faces
    for (label fI = 0; fI < nFaces; fI++)
    {
        const vector& n(F.normals()[fI]);
        for (label nG = 0; nG < nGaussPerFace; nG++)
        {
            scalar UL = U.minusValueOnFace(fI, nG);
            scalar UR = U.plusValueOnFace(fI, nG);

            vector FMinus = F.minusValueOnFace(fI, nG);
            vector FPlus = F.plusValueOnFace(fI, nG);

            calcIntermediateState
            (
                cellID,
                fI,
                nG,
                n,

                // output
                SL,
                SR,
                SStar,
                CL,
                CR
            );

            // Compute U*L and U*R based on equationType
            scalar UStarL = 0.0;
            scalar UStarR = 0.0;

            switch (eqnType_)
            {
                case dgFluxSolver::equationType::massTransport:
                {
                    UStarL = CL;
                    UStarR = CR;
                    break;
                }
                case dgFluxSolver::equationType::energyTransport:
                {
                    // Left and right states
                    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
                    const scalar rhoR = rhoF.plusValueOnFace(fI, nG);
                    const scalar rhoL = rhoF.minusValueOnFace(fI, nG);

                    const faceGaussField<scalar>& pF = p_.gaussFields()[cellID].faceField();
                    const scalar pR = pF.plusValueOnFace(fI, nG);
                    const scalar pL = pF.minusValueOnFace(fI, nG);

                    const faceGaussField<vector>& UF = U_.gaussFields()[cellID].faceField();
                    const vector URv = UF.plusValueOnFace(fI, nG);
                    const vector ULv = UF.minusValueOnFace(fI, nG);

                    // Normal velocities
                    const scalar UnR = (URv & n);
                    const scalar UnL = (ULv & n);

                    const scalar EL = UL; // assuming UL = rho*e + 0.5*rho*U^2
                    const scalar ER = UR;
                    UStarL = CL * (EL/rhoL + (SStar - UnL)*(SStar + pL/(rhoL*(SL - UnL + VSMALL))));
                    UStarR = CR * (ER/rhoR + (SStar - UnR)*(SStar + pR/(rhoR*(SR - UnR + VSMALL))));
                    break;
                }
                case dgFluxSolver::equationType::scalarTransport:
                {
                    // Left and right states
                    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
                    const scalar rhoR = rhoF.plusValueOnFace(fI, nG);
                    const scalar rhoL = rhoF.minusValueOnFace(fI, nG);

                    const scalar qL = UL / rhoL;
                    const scalar qR = UR / rhoR;
                    UStarL = CL * qL;
                    UStarR = CR * qR;
                    break;
                }
                default:
                    FatalErrorInFunction << "Unsupported equationType in HLLC." << nl << exit(FatalError);
            }

            // Compute FStarL and FStarR
            const scalar fL = (FMinus & n);
            const scalar fR = (FPlus & n);


            const scalar FStarL = fL + SL * (UStarL - UL);
            const scalar FStarR = fR + SR * (UStarR - UR);


            // HLLC logic (scalar flux fn)
            scalar fn = 0.0;
            if (0 <= SL)
                fn = fL;
            else if (SL <= 0 && 0 <= SStar)
                fn = FStarL;
            else if (SStar <= 0 && 0 <= SR)
                fn = FStarR;
            else // 0 >= SR
                fn = fR;

            // Back to Cartesian
            F.fluxOnFace(fI, nG) = fn * n;
        }
    }
}

void Foam::dgHLLCFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>& U
)
{
    scalar SL(0.0);
    scalar SR(0.0);
    scalar SStar(0.0);
    scalar CL(0.0);
    scalar CR(0.0);

    const label nFaces(F.nFaces());
    const label nGaussPerFace(F.nGaussPerFace());

    // Loop over faces
    for (label fI = 0; fI < nFaces; fI++)
    {
        const vector& n(F.normals()[fI]);
        for (label nG = 0; nG < nGaussPerFace; nG++)
        {
            vector UL = U.minusValueOnFace(fI, nG);
            vector UR = U.plusValueOnFace(fI, nG);

            tensor FMinus = F.minusValueOnFace(fI, nG);
            tensor FPlus = F.plusValueOnFace(fI, nG);

            calcIntermediateState
            (
                cellID,
                fI,
                nG,
                n,

                // output
                SL,
                SR,
                SStar,
                CL,
                CR
            );

            // Compute U*L and U*R based on equationType
            vector UStarL(Zero);
            vector UStarR(Zero);

            switch (eqnType_)
            {
                case dgFluxSolver::equationType::momentumTransport:
                {
                    const faceGaussField<vector>& UF   = U_.gaussFields()[cellID].faceField();

                    // Left (-) and right (+) states
                    const vector URv = UF.plusValueOnFace(fI, nG);
                    const vector ULv = UF.minusValueOnFace(fI, nG);

                    // normal and tangential component of velocity
                    vector t1, t2;
                    makeONB(n, t1, t2);

                    // decompose velocity
                    vector URvn, URvt1, URvt2, ULvn, ULvt1, ULvt2;
                    decomposeU(URv, n, URvn, URvt1, URvt2);
                    decomposeU(ULv, n, ULvn, ULvt1, ULvt2);

                    // reconstruct UStar
                    UStarL = CL*(SStar*n + ULvt1 + ULvt2);
                    UStarR = CR*(SStar*n + URvt1 + URvt2);
                    break;
                }
                default:
                    FatalErrorInFunction << "Unsupported equationType in HLLC." << nl << exit(FatalError);
            }

            // Compute FStarL and FStarR
            const vector fL = (FMinus & n);
            const vector fR = (FPlus & n);


            const vector FStarL = fL + SL * (UStarL - UL);
            const vector FStarR = fR + SR * (UStarR - UR);

            if (0 <= SL)
                F.fluxOnFace(fI, nG) = fL;
            else if (SL <= 0 && 0 <= SStar)
                F.fluxOnFace(fI, nG) = FStarL;
            else if (SStar <= 0 && 0 <= SR)
                F.fluxOnFace(fI, nG) = FStarR;
            else // 0 >= SR
                F.fluxOnFace(fI, nG) = fR;
        }
    }
}

void Foam::dgHLLCFluxSolver::reset()
{
    // Reset isStateComputed_ flag
    const label nFaces = mesh_.nFaces();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        isStateComputed_[fI] = false;
    }
}

} // End namespace Foam

// ************************************************************************* //

