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

#include "dgRoeFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

// #include "faceGaussField.H"
// #include "fieldsContext.H"

namespace Foam
{

defineTypeNameAndDebug(dgRoeFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgRoeFluxSolver, dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
Foam::dgRoeFluxSolver::dgRoeFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgThermo& thermo
)
:
    dgFluxSolver(name, dict, thermo)
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dgRoeFluxSolver::read(const dictionary& dict)
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

void Foam::dgHLLECFluxSolver::computeFlux
(
    const label gaussID,
    const vector& FPlus,   // N-side physical flux (Cartesian)
    const vector& FMinus,  // P-side physical flux (Cartesian)
    const scalar UR,       // conserved on N-side
    const scalar UL,       // conserved on P-side
    const vector& n,       // unit normal P -> N
    vector& flux           // output Cartesian numerical flux
) const
{
    if (!ctxPtr_)
    {
        FatalErrorInFunction
            << "fieldsContext is null. Call setContext(...) before computeFlux."
            << nl << exit(FatalError);
    }


    // 1) Lookup primitive fields
    const auto& rhoF = ctxPtr_->lookupFaceField<scalar>("rho");
    const auto& UF = ctxPtr_->lookupFaceField<vector>("U");
    const auto& pF = ctxPtr_->lookupFaceField<scalar>("p");


    const scalar rhoR = rhoF.plusValue(gaussID);
    const scalar rhoL = rhoF.minusValue(gaussID);
    const vector URv = UF.plusValue(gaussID);
    const vector ULv = UF.minusValue(gaussID);
    const scalar pR = pF.plusValue(gaussID);
    const scalar pL = pF.minusValue(gaussID);


    const scalar UnR = (URv & n);
    const scalar UnL = (ULv & n);


    const scalar TR = thermo_.eqnOfState().T(rhoR, pR);
    const scalar TL = thermo_.eqnOfState().T(rhoL, pL);
    const scalar aR = thermo_.thermo().a(TR);
    const scalar aL = thermo_.thermo().a(TL);


    // 2) Compute SL and SR
    scalar SL(0), SR(0);

    if (speedEst_ == seDavis)
    {
        // Davis/Einfeldt (robust/positivity-friendly)
        SL = min(UnL - aL, UnR - aR);
        SR = max(UnL + aL, UnR + aR);
    }
    else // seRoeEinfeldt
    {
        // Roe averages for U, H, a (less diffusive)
        const scalar sL = ::sqrt(max(rhoL, SMALL));
        const scalar sR = ::sqrt(max(rhoR, SMALL));
        const scalar denom = sL + sR + SMALL;

        const vector Uroe = (sL*ULv + sR*URv)/denom;
        const scalar UnRoe = (Uroe & n);

        const scalar Cp = thermo_.Cp();
        const scalar hL = Cp*TL;
        const scalar hR = Cp*TR;
        const scalar HL = hL + 0.5*magSqr(ULv);
        const scalar HR = hR + 0.5*magSqr(URv);
        const scalar HRoe = (sL*HL + sR*HR)/denom;

        const scalar gamma = thermo_.gamma();
        const scalar aRoe2 = max((gamma - 1.0)*(HRoe - 0.5*magSqr(Uroe)), 0.0);
        const scalar aRoe  = ::sqrt(aRoe2);

        SL = min(UnL - aL, UnRoe - aRoe);
        SR = max(UnR + aR, UnRoe + aRoe);
    }


    // 3) Compute S*
    const scalar SStar = 
        (pR - pL + rhoL * UnL * (SL - UnL) - rhoR * UnR * (SR - UnR))/
        (rhoL * (SL - UnL) - rhoR * (SR - UnR) + VSMALL);


    // 4) Compute Ck = rhoK*(Sk - uK)/(Sk - S*)
    const scalar CL = rhoL * (SL - UnL) / (SL - SStar + VSMALL);
    const scalar CR = rhoR * (SR - UnR) / (SR - SStar + VSMALL);


    // 5) Compute U*L and U*R based on equationType
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
        case dgFluxSolver::equationType::momentumTransport:
        {
            UStarL = CL * SStar;
            UStarR = CR * SStar;
            break;
        }
        case dgFluxSolver::equationType::energyTransport:
        {
            const scalar EL = UL; // assuming UL = rho*e + 0.5*rho*U^2
            const scalar ER = UR;
            UStarL = CL * (EL/rhoL + (SStar - UnL)*(SStar + pL/(rhoL*(SL - UnL + VSMALL))));
            UStarR = CR * (ER/rhoR + (SStar - UnR)*(SStar + pR/(rhoR*(SR - UnR + VSMALL))));
            break;
        }
        case dgFluxSolver::equationType::scalarTransport:
        {
            const scalar qL = UL / rhoL;
            const scalar qR = UR / rhoR;
            UStarL = CL * qL;
            UStarR = CR * qR;
            break;
        }
        default:
            FatalErrorInFunction << "Unsupported equationType in HLLC." << nl << exit(FatalError);
    }


    // 6) Compute FStarL and FStarR
    const scalar fL = (FMinus & n);
    const scalar fR = (FPlus & n);


    const scalar FStarL = fL + SL * (UStarL - UL);
    const scalar FStarR = fR + SR * (UStarR - UR);


    // 7) HLLC logic (scalar flux fn)
    scalar fn = 0.0;
    if (0 <= SL)
        fn = fL;
    else if (SL <= 0 && 0 <= SStar)
        fn = FStarL;
    else if (SStar <= 0 && 0 <= SR)
        fn = FStarR;
    else // 0 >= SR
        fn = fR;


    // 8) Back to Cartesian
    flux = fn * n;
}

} // End namespace Foam

// ************************************************************************* //

