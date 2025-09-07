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

// #include "faceGaussField.H"
// #include "fieldsContext.H"

namespace Foam
{

defineTypeNameAndDebug(dgHLLEFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgHLLEFluxSolver, dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
Foam::dgHLLEFluxSolver::dgHLLEFluxSolver
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

void Foam::dgHLLEFluxSolver::read(const dictionary& dict)
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

void Foam::dgHLLEFluxSolver::computeFlux
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
    // 0) Context check
    if (!ctxPtr_)
    {
        FatalErrorInFunction
            << "fieldsContext is null. Call setContext(...) before computeFlux."
            << nl << exit(FatalError);
    }

    // 1) Primitive fields at Gauss point
    const auto& rhoF = ctxPtr_->lookupFaceField<scalar>("rho");
    const auto& UF   = ctxPtr_->lookupFaceField<vector>("U");
    const auto& pF   = ctxPtr_->lookupFaceField<scalar>("p");

    const scalar rhoR = rhoF.plusValue(gaussID);   // N-side (Plus)
    const scalar rhoL = rhoF.minusValue(gaussID);  // P-side (Minus)
    const vector URv  = UF.plusValue(gaussID);
    const vector ULv  = UF.minusValue(gaussID);
    const scalar pR   = pF.plusValue(gaussID);
    const scalar pL   = pF.minusValue(gaussID);

    // 2) Normal velocities & sound speeds
    const scalar UnR = (URv & n);
    const scalar UnL = (ULv & n);

    const scalar TR = thermo_.eqnOfState().T(rhoR, pR);
    const scalar TL = thermo_.eqnOfState().T(rhoL, pL);
    const scalar aR = thermo_.thermo().a(TR);
    const scalar aL = thermo_.thermo().a(TL);

    // 3) Wave-speed bounds
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

    // 4) Physical normal fluxes (project F±)
    const scalar fR = (FPlus  & n);  // right  (N)
    const scalar fL = (FMinus & n);  // left   (P)

    // 5) HLLE normal flux with jump (UR - UL)
    scalar fn(0.0);
    if (SL >= 0.0)
    {
        fn = fL;
    }
    else if (SR <= 0.0)
    {
        fn = fR;
    }
    else
    {
        fn = (SR*fL - SL*fR + SL*SR*(UR - UL)) / (SR - SL + SMALL);
    }

    // 6) Back to Cartesian (tangential components vanish in 1D Riemann)
    flux = fn * n;
}

} // End namespace Foam

// ************************************************************************* //

