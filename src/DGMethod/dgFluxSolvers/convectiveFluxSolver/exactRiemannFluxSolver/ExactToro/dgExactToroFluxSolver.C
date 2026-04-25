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

#include "dgExactToroFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "dgThermoConservative.H"
#include "error.H"
#include <cmath>

namespace Foam
{
namespace
{

static constexpr scalar exactToroTol = 1e-6;
static constexpr label exactToroMaxIter = 20;

inline void calcRusanovFlux
(
    const scalar rhoL,
    const scalar rhoR,
    const vector& ULv,
    const vector& URv,
    const scalar EL,
    const scalar ER,
    const scalar pL,
    const scalar pR,
    const scalar gamma,
    const vector& n,
    scalar& massFlux,
    vector& momentumFlux,
    scalar& energyFlux
)
{
    const scalar rhoLSafe = max(rhoL, SMALL);
    const scalar rhoRSafe = max(rhoR, SMALL);
    const scalar pLSafe = max(pL, SMALL);
    const scalar pRSafe = max(pR, SMALL);
    const scalar gammaSafe = max(gamma, scalar(1.0 + SMALL));

    const scalar unL = (ULv & n);
    const scalar unR = (URv & n);
    const scalar aL = sqrt(max(gammaSafe*pLSafe/rhoLSafe, SMALL));
    const scalar aR = sqrt(max(gammaSafe*pRSafe/rhoRSafe, SMALL));
    const scalar sMax = max(mag(unL) + aL, mag(unR) + aR);

    const scalar massFluxL = rhoL*unL;
    const scalar massFluxR = rhoR*unR;

    const vector momentumFluxL = massFluxL*ULv + pL*n;
    const vector momentumFluxR = massFluxR*URv + pR*n;

    const scalar energyFluxL = unL*(EL + pL);
    const scalar energyFluxR = unR*(ER + pR);

    massFlux =
        0.5*(massFluxL + massFluxR)
      - 0.5*sMax*(rhoR - rhoL);

    momentumFlux =
        0.5*(momentumFluxL + momentumFluxR)
      - 0.5*sMax*(rhoR*URv - rhoL*ULv);

    energyFlux =
        0.5*(energyFluxL + energyFluxR)
      - 0.5*sMax*(ER - EL);
}

struct ToroPrimitiveState
{
    scalar rho;
    scalar u;
    scalar v;
    scalar w;
    scalar p;
};

inline scalar guessp
(
    const scalar gamma[9],
    const scalar rhoL,
    const scalar uL,
    const scalar pL,
    const scalar cL,
    const scalar rhoR,
    const scalar uR,
    const scalar pR,
    const scalar cR
)
{
    const scalar quser = 2.0;

    // Start from Toro's primitive-variable pressure estimate (PVRS).
    scalar cup = 0.25*(rhoL + rhoR)*(cL + cR);
    scalar ppv = 0.5*(pL + pR) + 0.5*(uL - uR)*cup;
    ppv = max(ppv, scalar(0));

    const scalar pmin = min(pL, pR);
    const scalar pmax = max(pL, pR);
    const scalar qmax = pmax/max(pmin, VSMALL);

    // If the states are not strongly separated, the PVRS estimate is enough.
    if (qmax <= quser && pmin <= ppv && ppv <= pmax)
    {
        return ppv;
    }

    // Rarefaction on both sides: use Toro's two-rarefaction estimate.
    if (ppv < pmin)
    {
        const scalar pq = pow(pL/pR, gamma[1]);
        const scalar um =
        (
            pq*uL/cL + uR/cR + gamma[4]*(pq - 1.0)
        ) / (pq/cL + 1.0/cR);

        const scalar ptL = 1.0 + gamma[7]*(uL - um)/cL;
        const scalar ptR = 1.0 + gamma[7]*(um - uR)/cR;

        return 0.5*
        (
            pL*pow(ptL, gamma[3]) + pR*pow(ptR, gamma[3])
        );
    }

    // Otherwise at least one side is a shock, so use the two-shock estimate.
    const scalar geL = sqrt((gamma[5]/rhoL)/(gamma[6]*pL + ppv));
    const scalar geR = sqrt((gamma[5]/rhoR)/(gamma[6]*pR + ppv));

    return (geL*pL + geR*pR - (uR - uL))/(geL + geR);
}


inline void prefun
(
    const scalar gamma[9],
    const scalar p,
    const scalar rho,
    const scalar pk,
    const scalar ck,
    scalar& f,
    scalar& fd
)
{
    if (p <= pk)
    {
        // Rarefaction branch of Toro's pressure function f(p).
        const scalar prat = p/pk;
        f = gamma[4]*ck*(pow(prat, gamma[1]) - 1.0);
        fd = pow(prat, -gamma[2])/(rho*ck);
        return;
    }

    // Shock branch of Toro's pressure function f(p).
    const scalar ak = gamma[5]/rho;
    const scalar bk = gamma[6]*pk;
    const scalar qrt = sqrt(ak/(bk + p));

    f = (p - pk)*qrt;
    fd = (1.0 - 0.5*(p - pk)/(bk + p))*qrt;
}


inline ToroPrimitiveState sampleToroState
(
    const scalar gamma,
    const scalar gammaArr[9],
    const scalar rhoL,
    const scalar uL,
    const scalar vL,
    const scalar wL,
    const scalar pL,
    const scalar cL,
    const scalar rhoR,
    const scalar uR,
    const scalar vR,
    const scalar wR,
    const scalar pR,
    const scalar cR,
    const scalar pStar,
    const scalar uStar
)
{
    const scalar S = 0.0;
    ToroPrimitiveState out{};

    // Sample the exact self-similar solution at x/t = 0, i.e. the interface.
    if (S <= uStar)
    {
        // The interface lies on the left side of the contact discontinuity.
        if (pStar <= pL)
        {
            // Left-moving rarefaction.
            const scalar shL = uL - cL;

            if (S <= shL)
            {
                // Sampling point is still in the original left state.
                out = {rhoL, uL, vL, wL, pL};
            }
            else
            {
                const scalar cmL = cL*pow(pStar/pL, gammaArr[1]);
                const scalar stL = uStar - cmL;

                if (S > stL)
                {
                    // Sampling point is inside the left star state.
                    out = {rhoL*pow(pStar/pL, 1.0/gamma), uStar, vL, wL, pStar};
                }
                else
                {
                    // Sampling point lies inside the left rarefaction fan.
                    const scalar c = gammaArr[5]*(cL + gammaArr[7]*(uL - S));
                    out =
                    {
                        rhoL*pow(c/cL, gammaArr[4]),
                        gammaArr[5]*(cL + gammaArr[7]*uL + S),
                        vL,
                        wL,
                        pL*pow(c/cL, gammaArr[3])
                    };
                }
            }
        }
        else
        {
            // Left-moving shock.
            const scalar pmL = pStar/pL;
            const scalar SL = uL - cL*sqrt(gammaArr[2]*pmL + gammaArr[1]);

            if (S <= SL)
            {
                // Sampling point is ahead of the left shock.
                out = {rhoL, uL, vL, wL, pL};
            }
            else
            {
                // Sampling point is behind the left shock, in the star region.
                out =
                {
                    rhoL*(pmL + gammaArr[6])/(pmL*gammaArr[6] + 1.0),
                    uStar,
                    vL,
                    wL,
                    pStar
                };
            }
        }
    }
    else
    {
        // The interface lies on the right side of the contact discontinuity.
        if (pStar > pR)
        {
            // Right-moving shock.
            const scalar pmR = pStar/pR;
            const scalar SR = uR + cR*sqrt(gammaArr[2]*pmR + gammaArr[1]);

            if (S >= SR)
            {
                // Sampling point is still in the original right state.
                out = {rhoR, uR, vR, wR, pR};
            }
            else
            {
                // Sampling point is in the right star state behind the shock.
                out =
                {
                    rhoR*(pmR + gammaArr[6])/(pmR*gammaArr[6] + 1.0),
                    uStar,
                    vR,
                    wR,
                    pStar
                };
            }
        }
        else
        {
            // Right-moving rarefaction.
            const scalar shR = uR + cR;

            if (S >= shR)
            {
                // Sampling point is ahead of the right rarefaction.
                out = {rhoR, uR, vR, wR, pR};
            }
            else
            {
                const scalar cmR = cR*pow(pStar/pR, gammaArr[1]);
                const scalar stR = uStar + cmR;

                if (S <= stR)
                {
                    // Sampling point is inside the right star state.
                    out = {rhoR*pow(pStar/pR, 1.0/gamma), uStar, vR, wR, pStar};
                }
                else
                {
                    // Sampling point lies inside the right rarefaction fan.
                    const scalar c = gammaArr[5]*(cR - gammaArr[7]*(uR - S));
                    out =
                    {
                        rhoR*pow(c/cR, gammaArr[4]),
                        gammaArr[5]*(-cR + gammaArr[7]*uR + S),
                        vR,
                        wR,
                        pR*pow(c/cR, gammaArr[3])
                    };
                }
            }
        }
    }

    return out;
}

} // End anonymous namespace

defineTypeNameAndDebug(dgExactToroFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgExactToroFluxSolver, dictionary);

dgExactToroFluxSolver::dgExactToroFluxSolver
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgFluxSolver(name, dict, mesh),
    thermo_
    (
        mesh.getFvMesh().lookupObject<dgThermoConservative>("dgThermoConservative")
    ),
    rho_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("rho")
    ),
    E_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("E")
    ),
    U_
    (
        mesh.getFvMesh().lookupObject<dgField<vector>>("U")
    ),
    p_
    (
        mesh.getFvMesh().lookupObject<dgField<scalar>>("p")
    ),
    warnedVariableGamma_(false),
    warnedVacuumFallback_(false),
    warnedNewtonFallback_(false)
{
    fType_ = dgFluxSolver::fluxType::convective;
    read(dict);
}


void dgExactToroFluxSolver::read(const dictionary&)
{}


void dgExactToroFluxSolver::calcGamma
(
    const label cellID,
    GaussField<scalar>& gamma
) const
{
    const GaussField<scalar>& TG = thermo_.T().gaussFields()[cellID];

    tmp<GaussField<scalar>> tCp = GaussField<scalar>::New(cellID, &mesh_);
    tmp<GaussField<scalar>> tCv = GaussField<scalar>::New(cellID, &mesh_);

    // Recover Cp, Cv, and gamma at the cell Gauss points from the thermo model.
    thermo_.thermo().calcCp(cellID, TG, tCp.ref());
    thermo_.thermo().calcCv(cellID, TG, tCv.ref());
    thermo_.thermo().calcGamma(cellID, tCp(), tCv(), gamma);
}


scalar dgExactToroFluxSolver::selectInterfaceGamma
(
    const scalar gammaL,
    const scalar gammaR
) const
{
    if (!warnedVariableGamma_ && mag(gammaL - gammaR) > 1e-10)
    {
        warnedVariableGamma_ = true;

        // Toro's exact derivation assumes a single constant gamma on both sides.
        WarningInFunction
            << "ExactToro assumes a calorically perfect gas with constant gamma. "
            << "Detected gammaL=" << gammaL
            << " and gammaR=" << gammaR
            << ". Using their arithmetic mean at the interface." << nl;
    }

    return 0.5*(gammaL + gammaR);
}


void dgExactToroFluxSolver::calcExactFlux
(
    const scalar rhoL,
    const scalar rhoR,
    const vector& ULv,
    const vector& URv,
    const scalar EL,
    const scalar ER,
    const scalar pL,
    const scalar pR,
    const scalar gammaL,
    const scalar gammaR,
    const vector& n,
    scalar& massFlux,
    vector& momentumFlux,
    scalar& energyFlux
) const
{
    // Protect the exact solver from unphysical states entering square-roots
    // or divisions during the 1-D Riemann solve.
    const scalar rhoLSafe = max(rhoL, SMALL);
    const scalar rhoRSafe = max(rhoR, SMALL);
    const scalar pLSafe = max(pL, SMALL);
    const scalar pRSafe = max(pR, SMALL);
    const scalar gamma = max(selectInterfaceGamma(gammaL, gammaR), scalar(1.0 + SMALL));

    // Build a local orthonormal basis so the exact solver only sees the
    // normal component, while tangential velocities are carried through.
    vector t1(Zero);
    vector t2(Zero);
    makeONB(n, t1, t2);

    const scalar uL = (ULv & n);
    const scalar vL = (ULv & t1);
    const scalar wL = (ULv & t2);

    const scalar uR = (URv & n);
    const scalar vR = (URv & t1);
    const scalar wR = (URv & t2);

    const scalar gammaArr[9] =
    {
        gamma,
        (gamma - 1.0)/(2.0*gamma),
        (gamma + 1.0)/(2.0*gamma),
        2.0*gamma/(gamma - 1.0),
        2.0/(gamma - 1.0),
        2.0/(gamma + 1.0),
        (gamma - 1.0)/(gamma + 1.0),
        0.5*(gamma - 1.0),
        gamma - 1.0
    };

    // Compute left/right acoustic speeds in the face-normal Riemann problem.
    const scalar cL = sqrt(max(gamma*pLSafe/rhoLSafe, SMALL));
    const scalar cR = sqrt(max(gamma*pRSafe/rhoRSafe, SMALL));

    // Toro's vacuum criterion: stop early if the input states would create vacuum.
    if (gammaArr[4]*(cL + cR) <= (uR - uL))
    {
        if (!warnedVacuumFallback_)
        {
            warnedVacuumFallback_ = true;

            WarningInFunction
                << "ExactToro detected an interface state satisfying Toro's "
                << "vacuum criterion and will fall back to local Rusanov flux."
                << nl
                << "rhoL=" << rhoL << ", rhoR=" << rhoR
                << ", pL=" << pL << ", pR=" << pR
                << ", unL=" << uL << ", unR=" << uR
                << ", gamma=" << gamma << nl;
        }

        calcRusanovFlux
        (
            rhoL, rhoR,
            ULv, URv,
            EL, ER,
            pL, pR,
            gamma,
            n,
            massFlux,
            momentumFlux,
            energyFlux
        );

        return;
    }

    // Use Toro's pressure guess as the starting point for Newton iteration.
    scalar pOld = guessp
    (
        gammaArr,
        rhoLSafe, uL, pLSafe, cL,
        rhoRSafe, uR, pRSafe, cR
    );

    const scalar uDiff = uR - uL;

    scalar pStar = pOld;
    scalar fL = 0.0;
    scalar fR = 0.0;
    scalar fLd = 0.0;
    scalar fRd = 0.0;
    scalar change = 0.0;
    label iter = 0;

    for (; iter < exactToroMaxIter; ++iter)
    {
        // Evaluate the exact left/right pressure functions and their derivatives.
        prefun(gammaArr, pOld, rhoLSafe, pLSafe, cL, fL, fLd);
        prefun(gammaArr, pOld, rhoRSafe, pRSafe, cR, fR, fRd);

        // Newton-Raphson update for the star-region pressure p*.
        pStar = pOld - (fL + fR + uDiff)/(fLd + fRd + VSMALL);
        pStar = max(pStar, scalar(exactToroTol));

        // Convergence monitor recommended by Toro.
        change = 2.0*mag((pStar - pOld)/(pStar + pOld + VSMALL));

        if (change <= exactToroTol)
        {
            break;
        }

        pOld = pStar;
    }

    if (iter >= exactToroMaxIter)
    {
        if (!warnedNewtonFallback_)
        {
            warnedNewtonFallback_ = true;

            WarningInFunction
                << "Newton-Raphson iteration failed to converge in ExactToro. "
                << "Falling back to local Rusanov flux." << nl
                << "rhoL=" << rhoL << ", rhoR=" << rhoR
                << ", pL=" << pL << ", pR=" << pR
                << ", unL=" << uL << ", unR=" << uR
                << ", gamma=" << gamma << nl;
        }

        calcRusanovFlux
        (
            rhoL, rhoR,
            ULv, URv,
            EL, ER,
            pL, pR,
            gamma,
            n,
            massFlux,
            momentumFlux,
            energyFlux
        );

        return;
    }

    // Once p* is known, the star-region normal velocity follows analytically.
    const scalar uStar = 0.5*(uL + uR + fR - fL);

    // Sample the exact solution at the interface to obtain the upwind state.
    const ToroPrimitiveState sampled = sampleToroState
    (
        gamma,
        gammaArr,
        rhoLSafe, uL, vL, wL, pLSafe, cL,
        rhoRSafe, uR, vR, wR, pRSafe, cR,
        pStar, uStar
    );

    const vector sampledU =
        sampled.u*n + sampled.v*t1 + sampled.w*t2;

    // Reconstruct total energy density from the sampled primitive state.
    const scalar sampledE =
        sampled.p/(gamma - 1.0)
      + 0.5*sampled.rho*
        (
            sqr(sampled.u) + sqr(sampled.v) + sqr(sampled.w)
        );

    // Convert the sampled star/interface state into Euler fluxes.
    massFlux = sampled.rho*sampled.u;
    momentumFlux = massFlux*sampledU + sampled.p*n;
    energyFlux = sampled.u*(sampledE + sampled.p);

    (void)EL;
    (void)ER;
}


void dgExactToroFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<vector>& F,
    const faceGaussField<scalar>&
)
{
    if
    (
        eqnType_ != dgFluxSolver::equationType::massTransport
     && eqnType_ != dgFluxSolver::equationType::energyTransport
     && eqnType_ != dgFluxSolver::equationType::scalarTransport
    )
    {
        FatalErrorInFunction
            << "Scalar ExactToro flux is only valid for mass or energy transport."
            << abort(FatalError);
    }

    if (eqnType_ == dgFluxSolver::equationType::scalarTransport)
    {
        FatalErrorInFunction
            << "ExactToro is an Euler Riemann solver and does not support "
            << "generic scalarTransport." << abort(FatalError);
    }

    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& EF = E_.gaussFields()[cellID].faceField();
    const faceGaussField<vector>& UF = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& pF = p_.gaussFields()[cellID].faceField();

    tmp<GaussField<scalar>> tGamma = GaussField<scalar>::New(cellID, &mesh_);
    calcGamma(cellID, tGamma.ref());
    const faceGaussField<scalar>& gammaF = tGamma().faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    // Solve one independent exact Riemann problem per face Gauss point.
    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            const scalar rhoMinus = rhoF.minusValueOnFace(fI, gI);
            const scalar rhoPlus = rhoF.plusValueOnFace(fI, gI);

            const vector UMinus = UF.minusValueOnFace(fI, gI);
            const vector UPlus = UF.plusValueOnFace(fI, gI);

            const scalar EMinus = EF.minusValueOnFace(fI, gI);
            const scalar EPlus = EF.plusValueOnFace(fI, gI);

            const scalar pMinus = pF.minusValueOnFace(fI, gI);
            const scalar pPlus = pF.plusValueOnFace(fI, gI);

            const scalar gammaMinus = gammaF.minusValueOnFace(fI, gI);
            const scalar gammaPlus = gammaF.plusValueOnFace(fI, gI);

            scalar massFlux = 0.0;
            vector momentumFlux(Zero);
            scalar energyFlux = 0.0;

            calcExactFlux
            (
                rhoMinus,
                rhoPlus,
                UMinus,
                UPlus,
                EMinus,
                EPlus,
                pMinus,
                pPlus,
                gammaMinus,
                gammaPlus,
                n,
                massFlux,
                momentumFlux,
                energyFlux
            );

            if (eqnType_ == dgFluxSolver::equationType::massTransport)
            {
                // Scalar mass equation stores rho*u_n projected on the face normal.
                F.fluxOnFace(fI, gI) = massFlux*n;
            }
            else
            {
                // Scalar energy equation stores (E + p)u_n projected on the face normal.
                F.fluxOnFace(fI, gI) = energyFlux*n;
            }
        }
    }
}


void dgExactToroFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>&
)
{
    if (eqnType_ != dgFluxSolver::equationType::momentumTransport)
    {
        FatalErrorInFunction
            << "Vector ExactToro flux is only valid for momentum transport."
            << abort(FatalError);
    }

    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
    const faceGaussField<vector>& UF = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& EF = E_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& pF = p_.gaussFields()[cellID].faceField();

    tmp<GaussField<scalar>> tGamma = GaussField<scalar>::New(cellID, &mesh_);
    calcGamma(cellID, tGamma.ref());
    const faceGaussField<scalar>& gammaF = tGamma().faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    // As above, solve an exact interface problem at every face quadrature point.
    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            const scalar rhoMinus = rhoF.minusValueOnFace(fI, gI);
            const scalar rhoPlus = rhoF.plusValueOnFace(fI, gI);

            const vector UMinus = UF.minusValueOnFace(fI, gI);
            const vector UPlus = UF.plusValueOnFace(fI, gI);

            const scalar EMinus = EF.minusValueOnFace(fI, gI);
            const scalar EPlus = EF.plusValueOnFace(fI, gI);

            const scalar pMinus = pF.minusValueOnFace(fI, gI);
            const scalar pPlus = pF.plusValueOnFace(fI, gI);

            const scalar gammaMinus = gammaF.minusValueOnFace(fI, gI);
            const scalar gammaPlus = gammaF.plusValueOnFace(fI, gI);

            scalar massFlux = 0.0;
            vector momentumFlux(Zero);
            scalar energyFlux = 0.0;

            calcExactFlux
            (
                rhoMinus,
                rhoPlus,
                UMinus,
                UPlus,
                EMinus,
                EPlus,
                pMinus,
                pPlus,
                gammaMinus,
                gammaPlus,
                n,
                massFlux,
                momentumFlux,
                energyFlux
            );

            // Momentum transport uses the full vector flux rho*u*u_n + p*n.
            F.fluxOnFace(fI, gI) = momentumFlux;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
