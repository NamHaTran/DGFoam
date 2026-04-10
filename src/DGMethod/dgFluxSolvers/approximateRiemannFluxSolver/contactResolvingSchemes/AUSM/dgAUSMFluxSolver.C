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

#include "dgAUSMFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "dgThermoConservative.H"
#include "error.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dgAUSMFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgAUSMFluxSolver, dictionary);


dgAUSMFluxSolver::dgAUSMFluxSolver
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
    variant_(variantType::ausm2),
    beta_(0.125),
    alpha_(0.1875),
    Kp_(0.25),
    Ku_(0.75),
    sigma_(1.0),
    Mco_(0.01),
    usePressureMachCorrection_(true),
    useVelocityPressureCorrection_(true)
{
    read(dict);
}


void dgAUSMFluxSolver::setVariant(const variantType variant)
{
    variant_ = variant;

    beta_ = 0.0;
    alpha_ = 0.0;
    Kp_ = 0.0;
    Ku_ = 0.0;
    sigma_ = 1.0;
    Mco_ = 0.01;
    usePressureMachCorrection_ = false;
    useVelocityPressureCorrection_ = false;

    switch (variant_)
    {
        case variantType::ausm0:
            break;

        case variantType::ausm1:
            beta_ = 0.125;
            alpha_ = 0.1875;
            break;

        case variantType::ausm2:
            beta_ = 0.125;
            alpha_ = 0.1875;
            Kp_ = 0.25;
            Ku_ = 0.75;
            sigma_ = 1.0;
            usePressureMachCorrection_ = true;
            useVelocityPressureCorrection_ = true;
            break;

        case variantType::ausm3:
            beta_ = 0.125;
            alpha_ = 0.1875;
            Kp_ = 0.25;
            Ku_ = 0.75;
            sigma_ = 1.0;
            Mco_ = 0.01;
            usePressureMachCorrection_ = true;
            useVelocityPressureCorrection_ = true;
            break;
    }
}


void dgAUSMFluxSolver::read(const dictionary& dict)
{
    const dictionary* readDict = &dict;

    if (dict.found("fluxSolversCoeffs"))
    {
        const dictionary& coeffsDict = dict.subDict("fluxSolversCoeffs");
        const word coeffKey("AUSMCoeffs");

        if (coeffsDict.found(coeffKey))
        {
            readDict = &coeffsDict.subDict(coeffKey);
        }
    }

    const word variantName = readDict->lookupOrDefault<word>("variant", "AUSM2");

    if (variantName == "AUSM0" || variantName == "ausm0")
    {
        setVariant(variantType::ausm0);
    }
    else if (variantName == "AUSM1" || variantName == "ausm1")
    {
        setVariant(variantType::ausm1);
    }
    else if (variantName == "AUSM2" || variantName == "ausm2")
    {
        setVariant(variantType::ausm2);
    }
    else if (variantName == "AUSM3" || variantName == "ausm3")
    {
        setVariant(variantType::ausm3);
    }
    else
    {
        WarningInFunction
            << "Unknown AUSM variant \"" << variantName
            << "\". Falling back to AUSM2." << nl;
        setVariant(variantType::ausm2);
    }

    if (variant_ == variantType::ausm3)
    {
        Mco_ = readDict->lookupOrDefault<scalar>("Mco", 0.01);
    }
}


scalar dgAUSMFluxSolver::M1Function
(
    const label side,
    const scalar M
) const
{
    if (side == 0)
    {
        return 0.5*(M + mag(M));
    }

    return 0.5*(M - mag(M));
}


scalar dgAUSMFluxSolver::M2Function
(
    const label side,
    const scalar M
) const
{
    if (side == 0)
    {
        return 0.25*sqr(M + 1.0);
    }

    return -0.25*sqr(M - 1.0);
}


scalar dgAUSMFluxSolver::M4Function
(
    const label side,
    const scalar beta,
    const scalar M
) const
{
    if (mag(M) >= 1.0)
    {
        return M1Function(side, M);
    }

    scalar out = M2Function(side, M);

    if (side == 0)
    {
        out *= 1.0 - 16.0*beta*M2Function(1, M);
    }
    else
    {
        out *= 1.0 + 16.0*beta*M2Function(0, M);
    }

    return out;
}


scalar dgAUSMFluxSolver::P5Function
(
    const label side,
    const scalar alpha,
    const scalar M
) const
{
    if (mag(M) >= 1.0)
    {
        return (1.0/(M + sign(M)*VSMALL))*M1Function(side, M);
    }

    scalar out = M2Function(side, M);

    if (side == 0)
    {
        out *= (2.0 - M) - 16.0*alpha*M*M2Function(1, M);
    }
    else
    {
        out *= (-2.0 - M) + 16.0*alpha*M*M2Function(0, M);
    }

    return out;
}


void dgAUSMFluxSolver::calcAUSMFlux
(
    const scalar rhoL,
    const scalar rhoR,
    const vector& ULv,
    const vector& URv,
    const scalar EL,
    const scalar ER,
    const scalar pL,
    const scalar pR,
    const scalar aL,
    const scalar aR,
    const vector& n,
    scalar& massFlux,
    vector& momentumFlux,
    scalar& energyFlux
) const
{
    const scalar rhoLSafe = max(rhoL, SMALL);
    const scalar rhoRSafe = max(rhoR, SMALL);
    const scalar pLSafe = max(pL, SMALL);
    const scalar pRSafe = max(pR, SMALL);
    const scalar aLSafe = max(aL, SMALL);
    const scalar aRSafe = max(aR, SMALL);

    vector t1(Zero);
    vector t2(Zero);
    makeONB(n, t1, t2);

    const scalar uL = (ULv & n);
    const scalar vL = (ULv & t1);
    const scalar wL = (ULv & t2);

    const scalar uR = (URv & n);
    const scalar vR = (URv & t1);
    const scalar wR = (URv & t2);

    const scalar cA = max(0.5*(aLSafe + aRSafe), scalar(SMALL));
    const scalar ML = uL/cA;
    const scalar MR = uR/cA;

    scalar Mbar = M4Function(0, beta_, ML) + M4Function(1, beta_, MR);

    scalar pbar =
        pLSafe*P5Function(0, alpha_, ML)
      + pRSafe*P5Function(1, alpha_, MR);

    if (usePressureMachCorrection_)
    {
        const scalar Mtilde = 0.5*(sqr(ML) + sqr(MR));
        const scalar rhoA = 0.5*(rhoLSafe + rhoRSafe);
        scalar KpEff = Kp_;

        if (variant_ == variantType::ausm3)
        {
            const scalar Mo =
                sqrt(min(scalar(1.0), max(Mtilde, sqr(Mco_))));
            const scalar fa = max(Mo*(2.0 - Mo), scalar(VSMALL));
            KpEff /= fa;
        }

        const scalar Mp =
            -KpEff*((pRSafe - pLSafe)/(rhoA*sqr(cA) + VSMALL))
           *max(1.0 - sigma_*Mtilde, scalar(0.0));

        Mbar += Mp;
    }

    if (useVelocityPressureCorrection_)
    {
        const scalar rhoA = 0.5*(rhoLSafe + rhoRSafe);

        const scalar pu =
            -2.0*Ku_*rhoA*sqr(cA)*(MR - ML)
           *P5Function(0, alpha_, ML)
           *P5Function(1, alpha_, MR);

        pbar += pu;
    }

    if (Mbar >= 0.0)
    {
        massFlux = cA*Mbar*rhoLSafe;
        momentumFlux =
            (cA*Mbar*rhoLSafe*uL + pbar)*n
          + (cA*Mbar*rhoLSafe*vL)*t1
          + (cA*Mbar*rhoLSafe*wL)*t2;
        energyFlux = cA*Mbar*(EL + pLSafe);
    }
    else
    {
        massFlux = cA*Mbar*rhoRSafe;
        momentumFlux =
            (cA*Mbar*rhoRSafe*uR + pbar)*n
          + (cA*Mbar*rhoRSafe*vR)*t1
          + (cA*Mbar*rhoRSafe*wR)*t2;
        energyFlux = cA*Mbar*(ER + pRSafe);
    }
}


void dgAUSMFluxSolver::computeFlux
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
            << "Scalar AUSM flux is only valid for mass or energy transport."
            << abort(FatalError);
    }

    if (eqnType_ == dgFluxSolver::equationType::scalarTransport)
    {
        FatalErrorInFunction
            << "AUSM is an Euler-system Riemann solver and does not support "
            << "generic scalarTransport." << abort(FatalError);
    }

    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& EF = E_.gaussFields()[cellID].faceField();
    const faceGaussField<vector>& UF = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& pF = p_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& aF = thermo_.a().gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            scalar massFlux = 0.0;
            vector momentumFlux(Zero);
            scalar energyFlux = 0.0;

            calcAUSMFlux
            (
                rhoF.minusValueOnFace(fI, gI),
                rhoF.plusValueOnFace(fI, gI),
                UF.minusValueOnFace(fI, gI),
                UF.plusValueOnFace(fI, gI),
                EF.minusValueOnFace(fI, gI),
                EF.plusValueOnFace(fI, gI),
                pF.minusValueOnFace(fI, gI),
                pF.plusValueOnFace(fI, gI),
                aF.minusValueOnFace(fI, gI),
                aF.plusValueOnFace(fI, gI),
                n,
                massFlux,
                momentumFlux,
                energyFlux
            );

            if (eqnType_ == dgFluxSolver::equationType::massTransport)
            {
                F.fluxOnFace(fI, gI) = massFlux*n;
            }
            else
            {
                F.fluxOnFace(fI, gI) = energyFlux*n;
            }
        }
    }
}


void dgAUSMFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>&
)
{
    if (eqnType_ != dgFluxSolver::equationType::momentumTransport)
    {
        FatalErrorInFunction
            << "Vector AUSM flux is only valid for momentum transport."
            << abort(FatalError);
    }

    const faceGaussField<scalar>& rhoF = rho_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& EF = E_.gaussFields()[cellID].faceField();
    const faceGaussField<vector>& UF = U_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& pF = p_.gaussFields()[cellID].faceField();
    const faceGaussField<scalar>& aF = thermo_.a().gaussFields()[cellID].faceField();

    const label nFaces = F.nFaces();
    const label nGauss = F.nGaussPerFace();

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            scalar massFlux = 0.0;
            vector momentumFlux(Zero);
            scalar energyFlux = 0.0;

            calcAUSMFlux
            (
                rhoF.minusValueOnFace(fI, gI),
                rhoF.plusValueOnFace(fI, gI),
                UF.minusValueOnFace(fI, gI),
                UF.plusValueOnFace(fI, gI),
                EF.minusValueOnFace(fI, gI),
                EF.plusValueOnFace(fI, gI),
                pF.minusValueOnFace(fI, gI),
                pF.plusValueOnFace(fI, gI),
                aF.minusValueOnFace(fI, gI),
                aF.plusValueOnFace(fI, gI),
                n,
                massFlux,
                momentumFlux,
                energyFlux
            );

            F.fluxOnFace(fI, gI) = momentumFlux;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
