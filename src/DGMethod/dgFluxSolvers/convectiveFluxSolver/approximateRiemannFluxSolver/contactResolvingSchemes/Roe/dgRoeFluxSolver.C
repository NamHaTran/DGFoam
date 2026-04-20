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

#include "dgRoeFluxSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "dgThermoConservative.H"
#include "error.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dgRoeFluxSolver, 0);
addToRunTimeSelectionTable(dgFluxSolver, dgRoeFluxSolver, dictionary);


dgRoeFluxSolver::dgRoeFluxSolver
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
    warnedVariableGamma_(false)
{
    fType_ = dgFluxSolver::fluxType::convective;
    read(dict);
}


void dgRoeFluxSolver::read(const dictionary&)
{}


void dgRoeFluxSolver::calcGamma
(
    const label cellID,
    GaussField<scalar>& gamma
) const
{
    const GaussField<scalar>& TG = thermo_.T().gaussFields()[cellID];

    tmp<GaussField<scalar>> tCp = GaussField<scalar>::New(cellID, &mesh_);
    tmp<GaussField<scalar>> tCv = GaussField<scalar>::New(cellID, &mesh_);

    thermo_.thermo().calcCp(cellID, TG, tCp.ref());
    thermo_.thermo().calcCv(cellID, TG, tCv.ref());
    thermo_.thermo().calcGamma(cellID, tCp(), tCv(), gamma);
}


scalar dgRoeFluxSolver::selectInterfaceGamma
(
    const scalar gammaL,
    const scalar gammaR
) const
{
    if (!warnedVariableGamma_ && mag(gammaL - gammaR) > 1e-10)
    {
        warnedVariableGamma_ = true;

        WarningInFunction
            << "Roe uses a single interface gamma in the linearized Euler "
            << "system. Detected gammaL=" << gammaL
            << " and gammaR=" << gammaR
            << ". Using their arithmetic mean at the interface." << nl;
    }

    return 0.5*(gammaL + gammaR);
}


void dgRoeFluxSolver::calcRoeFlux
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
    const scalar rhoLSafe = max(rhoL, SMALL);
    const scalar rhoRSafe = max(rhoR, SMALL);
    const scalar pLSafe = max(pL, SMALL);
    const scalar pRSafe = max(pR, SMALL);
    const scalar gamma =
        max(selectInterfaceGamma(gammaL, gammaR), scalar(1.0 + SMALL));

    vector t1(Zero);
    vector t2(Zero);
    makeONB(n, t1, t2);

    const scalar uL = (ULv & n);
    const scalar vL = (ULv & t1);
    const scalar wL = (ULv & t2);

    const scalar uR = (URv & n);
    const scalar vR = (URv & t1);
    const scalar wR = (URv & t2);

    const scalar rhouL = rhoLSafe*uL;
    const scalar rhovL = rhoLSafe*vL;
    const scalar rhowL = rhoLSafe*wL;

    const scalar rhouR = rhoRSafe*uR;
    const scalar rhovR = rhoRSafe*vR;
    const scalar rhowR = rhoRSafe*wR;

    const scalar HL = (EL + pLSafe)/rhoLSafe;
    const scalar HR = (ER + pRSafe)/rhoRSafe;

    const scalar srL = sqrt(rhoLSafe);
    const scalar srR = sqrt(rhoRSafe);
    const scalar srLR = srL + srR + VSMALL;

    const scalar uRoe = (srL*uL + srR*uR)/srLR;
    const scalar vRoe = (srL*vL + srR*vR)/srLR;
    const scalar wRoe = (srL*wL + srR*wR)/srLR;
    const scalar hRoe = (srL*HL + srR*HR)/srLR;

    const scalar URoe2 = sqr(uRoe) + sqr(vRoe) + sqr(wRoe);
    const scalar cRoe2 =
        max((gamma - 1.0)*(hRoe - 0.5*URoe2), scalar(SMALL));
    const scalar cRoe = sqrt(cRoe2);

    const scalar k[5][5] =
    {
        {1.0, uRoe - cRoe, vRoe, wRoe, hRoe - uRoe*cRoe},
        {1.0, uRoe,        vRoe, wRoe, 0.5*URoe2},
        {0.0, 0.0,         1.0,  0.0,  vRoe},
        {0.0, 0.0,         0.0,  1.0,  wRoe},
        {1.0, uRoe + cRoe, vRoe, wRoe, hRoe + uRoe*cRoe}
    };

    const scalar jump[5] =
    {
        rhoRSafe - rhoLSafe,
        rhouR - rhouL,
        rhovR - rhovL,
        rhowR - rhowL,
        ER - EL
    };

    const scalar jumpBar =
        jump[4]
      - (jump[2] - vRoe*jump[0])*vRoe
      - (jump[3] - wRoe*jump[0])*wRoe;

    scalar alpha[5];

    alpha[1] =
        (gamma - 1.0)
       *(
            jump[0]*(hRoe - uRoe*uRoe)
          + uRoe*jump[1]
          - jumpBar
        )/(cRoe2 + VSMALL);

    alpha[0] =
        (
            jump[0]*(uRoe + cRoe)
          - jump[1]
          - cRoe*alpha[1]
        )/(2.0*cRoe + VSMALL);

    alpha[4] = jump[0] - (alpha[0] + alpha[1]);
    alpha[2] = jump[2] - vRoe*jump[0];
    alpha[3] = jump[3] - wRoe*jump[0];

    scalar rhof = 0.5*(rhoLSafe*uL + rhoRSafe*uR);
    scalar rhouf = 0.5*(pLSafe + rhoLSafe*uL*uL + pRSafe + rhoRSafe*uR*uR);
    scalar rhovf = 0.5*(rhoLSafe*uL*vL + rhoRSafe*uR*vR);
    scalar rhowf = 0.5*(rhoLSafe*uL*wL + rhoRSafe*uR*wR);
    scalar Ef = 0.5*(uL*(EL + pLSafe) + uR*(ER + pRSafe));

    const scalar lambda[5] =
    {
        mag(uRoe - cRoe),
        mag(uRoe),
        mag(uRoe),
        mag(uRoe),
        mag(uRoe + cRoe)
    };

    for (label i = 0; i < 5; ++i)
    {
        const scalar coeff = 0.5*alpha[i]*lambda[i];

        rhof -= coeff*k[i][0];
        rhouf -= coeff*k[i][1];
        rhovf -= coeff*k[i][2];
        rhowf -= coeff*k[i][3];
        Ef -= coeff*k[i][4];
    }

    massFlux = rhof;
    momentumFlux = rhouf*n + rhovf*t1 + rhowf*t2;
    energyFlux = Ef;
}


void dgRoeFluxSolver::computeFlux
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
            << "Scalar Roe flux is only valid for mass or energy transport."
            << abort(FatalError);
    }

    if (eqnType_ == dgFluxSolver::equationType::scalarTransport)
    {
        FatalErrorInFunction
            << "Roe is an Euler-system Riemann solver and does not support "
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

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            scalar massFlux = 0.0;
            vector momentumFlux(Zero);
            scalar energyFlux = 0.0;

            calcRoeFlux
            (
                rhoF.minusValueOnFace(fI, gI),
                rhoF.plusValueOnFace(fI, gI),
                UF.minusValueOnFace(fI, gI),
                UF.plusValueOnFace(fI, gI),
                EF.minusValueOnFace(fI, gI),
                EF.plusValueOnFace(fI, gI),
                pF.minusValueOnFace(fI, gI),
                pF.plusValueOnFace(fI, gI),
                gammaF.minusValueOnFace(fI, gI),
                gammaF.plusValueOnFace(fI, gI),
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


void dgRoeFluxSolver::computeFlux
(
    const label cellID,
    faceGaussField<tensor>& F,
    const faceGaussField<vector>&
)
{
    if (eqnType_ != dgFluxSolver::equationType::momentumTransport)
    {
        FatalErrorInFunction
            << "Vector Roe flux is only valid for momentum transport."
            << abort(FatalError);
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

    for (label fI = 0; fI < nFaces; ++fI)
    {
        const vector& n = F.normals()[fI];

        for (label gI = 0; gI < nGauss; ++gI)
        {
            scalar massFlux = 0.0;
            vector momentumFlux(Zero);
            scalar energyFlux = 0.0;

            calcRoeFlux
            (
                rhoF.minusValueOnFace(fI, gI),
                rhoF.plusValueOnFace(fI, gI),
                UF.minusValueOnFace(fI, gI),
                UF.plusValueOnFace(fI, gI),
                EF.minusValueOnFace(fI, gI),
                EF.plusValueOnFace(fI, gI),
                pF.minusValueOnFace(fI, gI),
                pF.plusValueOnFace(fI, gI),
                gammaF.minusValueOnFace(fI, gI),
                gammaF.plusValueOnFace(fI, gI),
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
