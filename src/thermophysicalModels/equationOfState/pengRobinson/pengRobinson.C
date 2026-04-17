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

#include "pengRobinson.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "error.H"

#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(pengRobinson, 0);
addToRunTimeSelectionTable(eqnOfState, pengRobinson, dictionary);

namespace
{

scalar readRequiredScalar
(
    const dictionary& dict,
    const word& primaryKey,
    const word& aliasKey = word::null
)
{
    if (dict.found(primaryKey))
    {
        return readScalar(dict.lookup(primaryKey));
    }

    if (!aliasKey.empty() && dict.found(aliasKey))
    {
        return readScalar(dict.lookup(aliasKey));
    }

    if (!aliasKey.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Required entry '" << primaryKey << "' (or alias '"
            << aliasKey << "') not found."
            << exit(FatalIOError);
    }

    FatalIOErrorInFunction(dict)
        << "Required entry '" << primaryKey << "' not found."
        << exit(FatalIOError);

    return Zero;
}

} // End anonymous namespace


pengRobinson::pengRobinson
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    eqnOfState(name, dict, mesh),
    R_(Zero),
    Tc_(Zero),
    Pc_(Zero),
    omega_(Zero),
    fw_(Zero),
    aCoeff_(Zero),
    bCoeff_(Zero),
    cpIdeal_(Zero),
    cvIdeal_(Zero),
    gammaIdeal_(Zero)
{
    if (!dict.found("specie"))
    {
        FatalIOErrorInFunction(dict)
            << "Peng-Robinson EOS requires a 'specie' sub-dictionary."
            << exit(FatalIOError);
    }

    if (!dict.found("thermodynamics"))
    {
        FatalIOErrorInFunction(dict)
            << "Peng-Robinson EOS requires a 'thermodynamics' sub-dictionary "
            << "with a reference ideal-gas Cp."
            << exit(FatalIOError);
    }

    if (!dict.found("equationOfState"))
    {
        FatalIOErrorInFunction(dict)
            << "Peng-Robinson EOS requires an 'equationOfState' "
            << "sub-dictionary."
            << exit(FatalIOError);
    }

    const dictionary& specieDict = dict.subDict("specie");
    const dictionary& thermoDict = dict.subDict("thermodynamics");
    const dictionary& eosDict = dict.subDict("equationOfState");

    const scalar molWeight = readRequiredScalar(specieDict, "molWeight");
    cpIdeal_ = readRequiredScalar(thermoDict, "Cp");
    Tc_ = readRequiredScalar(eosDict, "Tcrit", "Tc");
    Pc_ = readRequiredScalar(eosDict, "Pcrit", "Pc");
    omega_ = readRequiredScalar(eosDict, "AcentricFactor", "omega");

    if (molWeight <= SMALL)
    {
        FatalIOErrorInFunction(specieDict)
            << "Invalid 'molWeight' <= 0 for Peng-Robinson EOS."
            << exit(FatalIOError);
    }

    if (Tc_ <= SMALL || Pc_ <= SMALL)
    {
        FatalIOErrorInFunction(eosDict)
            << "Peng-Robinson requires positive critical temperature and "
            << "pressure."
            << exit(FatalIOError);
    }

    R_ = constant::thermodynamic::RR / molWeight;
    cvIdeal_ = cpIdeal_ - R_;

    if (cvIdeal_ <= SMALL)
    {
        FatalIOErrorInFunction(thermoDict)
            << "Computed Cv = Cp - R <= 0 for Peng-Robinson EOS. "
            << "Check 'Cp' and 'molWeight'."
            << exit(FatalIOError);
    }

    gammaIdeal_ = cpIdeal_ / cvIdeal_;
    fw_ = 0.37464 + 1.54226*omega_ - 0.26992*omega_*omega_;
    aCoeff_ = 0.45724*R_*R_*Tc_*Tc_/Pc_;
    bCoeff_ = 0.07780*R_*Tc_/Pc_;
}


scalar pengRobinson::alpha(const scalar T) const
{
    const scalar sqrtAlpha = 1.0 + fw_*(1.0 - std::sqrt(max(T, SMALL)/Tc_));
    return sqrtAlpha*sqrtAlpha;
}


scalar pengRobinson::logTerm(const scalar rho) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar oneOverRho = 1.0/rhoSafe;
    const scalar sqrtTwo = std::sqrt(2.0);
    const scalar numerator = oneOverRho + bCoeff_ - bCoeff_*sqrtTwo;
    const scalar denominator = oneOverRho + bCoeff_ + bCoeff_*sqrtTwo;

    return std::log(max(numerator/denominator, SMALL));
}


scalar pengRobinson::calcPFromRhoT
(
    const scalar rho,
    const scalar T
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar TSafe = max(T, SMALL);
    const scalar oneOverRho = 1.0/rhoSafe;
    const scalar denominator = oneOverRho*oneOverRho
        + 2.0*bCoeff_*oneOverRho
        - bCoeff_*bCoeff_;

    return
        R_*TSafe/(oneOverRho - bCoeff_)
      - aCoeff_*alpha(TSafe)/max(denominator, SMALL);
}


scalar pengRobinson::calcTemperatureFromRhoP
(
    const scalar rho,
    const scalar p
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar denom = 1.0/(rhoSafe*rhoSafe)
        + 2.0*bCoeff_/rhoSafe
        - bCoeff_*bCoeff_;

    const scalar A =
        R_/(1.0/rhoSafe - bCoeff_)
      - (aCoeff_*fw_*fw_)/(max(denom, SMALL)*Tc_);

    const scalar B =
        2.0*aCoeff_/max(denom, SMALL)
       *fw_*(1.0 + fw_)/std::sqrt(Tc_);

    const scalar C =
        -aCoeff_*(1.0 + fw_)*(1.0 + fw_)/max(denom, SMALL)
      - p;

    const scalar disc = max(B*B - 4.0*A*C, scalar(0));
    const scalar sqrtT = (-B + std::sqrt(disc))/(2.0*A + VSMALL);

    return max(sqrtT*sqrtT, SMALL);
}


scalar pengRobinson::calcTFromPRho
(
    const scalar p,
    const scalar rho
) const
{
    return calcTemperatureFromRhoP(rho, p);
}


scalar pengRobinson::calcRhoFromPT
(
    const scalar p,
    const scalar T
) const
{
    const scalar TSafe = max(T, SMALL);
    const scalar A =
        aCoeff_*alpha(TSafe)*p/(R_*R_*TSafe*TSafe);
    const scalar B = bCoeff_*p/(R_*TSafe);

    const scalar k1 = B - 1.0;
    const scalar k2 = A - 2.0*B - 3.0*B*B;
    const scalar k3 = -A*B + B*B + B*B*B;

    scalar Z = 1.0;
    scalar residual = GREAT;

    for (label iter = 0; iter < 100; ++iter)
    {
        const scalar f = Z*Z*Z + k1*Z*Z + k2*Z + k3;
        const scalar df = 3.0*Z*Z + 2.0*k1*Z + k2;

        if (mag(df) <= VSMALL)
        {
            break;
        }

        residual = f/df;
        Z -= residual;

        if (mag(residual) <= 1e-6)
        {
            break;
        }
    }

    if (mag(residual) > 1e-6)
    {
        WarningInFunction
            << "Peng-Robinson rho(p,T) iteration did not converge tightly; "
            << "using the latest compressibility estimate Z=" << Z << nl;
    }

    return p/(max(Z, SMALL)*R_*TSafe);
}


scalar pengRobinson::calcTFromRhoE
(
    const scalar rho,
    const scalar e
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar sqrtTwo = std::sqrt(2.0);
    const scalar logValue = logTerm(rhoSafe);
    const scalar A = cvIdeal_;
    const scalar factor1 = aCoeff_/(2.0*bCoeff_*sqrtTwo);
    const scalar factor2 = 1.0 + fw_;
    const scalar B =
        -factor1*logValue/std::sqrt(Tc_)*fw_*factor2;
    const scalar C = factor1*logValue*factor2*factor2 - e;
    const scalar disc = max(B*B - 4.0*A*C, scalar(0));
    const scalar sqrtT = (std::sqrt(disc) - B)/(2.0*A + VSMALL);

    return max(sqrtT*sqrtT, SMALL);
}


scalar pengRobinson::calcPFromRhoE
(
    const scalar rho,
    const scalar e
) const
{
    return calcPFromRhoT(rho, calcTFromRhoE(rho, e));
}


scalar pengRobinson::calcEFromRhoP
(
    const scalar rho,
    const scalar p
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar T = calcTemperatureFromRhoP(rhoSafe, p);
    const scalar sqrtAlpha = std::sqrt(alpha(T));
    const scalar sqrtTr = std::sqrt(max(T, SMALL)/Tc_);
    const scalar sqrtEight = std::sqrt(8.0);

    return
        cvIdeal_*T
      + aCoeff_/(bCoeff_*sqrtEight)*logTerm(rhoSafe)
       *(alpha(T) + sqrtAlpha*fw_*sqrtTr);
}


scalar pengRobinson::calcEFromRhoT
(
    const scalar rho,
    const scalar T
) const
{
    return calcEFromRhoP(rho, calcPFromRhoT(rho, T));
}


scalar pengRobinson::calcDPDeRho
(
    const scalar rho,
    const scalar e
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar T = calcTFromRhoE(rhoSafe, e);
    const scalar denom = 1.0/(rhoSafe*rhoSafe)
        + 2.0*bCoeff_/rhoSafe
        - bCoeff_*bCoeff_;
    const scalar sqrtAlpha = std::sqrt(alpha(T));
    const scalar dPdT =
        R_/(1.0/rhoSafe - bCoeff_)
      + aCoeff_/(std::sqrt(max(T*Tc_, SMALL)))
       *fw_*sqrtAlpha/max(denom, SMALL);

    const scalar cv =
        cvIdeal_
      - aCoeff_/(2.0*bCoeff_*std::sqrt(8.0))
       *logTerm(rhoSafe)
       *(fw_*(1.0 + fw_))/std::sqrt(max(T*Tc_, SMALL));

    return dPdT/max(cv, SMALL);
}


scalar pengRobinson::calcDPDrhoE
(
    const scalar rho,
    const scalar e
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar T = calcTFromRhoE(rhoSafe, e);
    const scalar dPde = calcDPDeRho(rhoSafe, e);
    const scalar denom = 1.0/(rhoSafe*rhoSafe)
        + 2.0*bCoeff_/rhoSafe
        - bCoeff_*bCoeff_;
    const scalar alphaT = alpha(T);
    const scalar dPdrhoT =
        R_*T/sqr(1.0 - bCoeff_*rhoSafe)
      - 2.0*aCoeff_*alphaT*rhoSafe*(1.0 + bCoeff_*rhoSafe)
       /sqr(max(denom*rhoSafe*rhoSafe, SMALL));

    const scalar dEdrhoT =
        -aCoeff_*std::sqrt(alphaT)*(1.0 + fw_)
       /max(denom*rhoSafe*rhoSafe, SMALL);

    return dPdrhoT - dPde*dEdrhoT;
}


scalar pengRobinson::calcAFromRhoE
(
    const scalar rho,
    const scalar e
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar p = calcPFromRhoE(rhoSafe, e);
    const scalar dpDe = calcDPDeRho(rhoSafe, e);
    const scalar dpDrho = calcDPDrhoE(rhoSafe, e);
    const scalar h = e + p/rhoSafe;
    const scalar chi = dpDrho - e/rhoSafe*dpDe;
    const scalar kappa = dpDe/rhoSafe;

    return std::sqrt(max(chi + kappa*h, SMALL));
}


vector pengRobinson::calcGradPFromRhoT
(
    const scalar rho,
    const scalar T,
    const vector& gradRho,
    const vector& gradT
) const
{
    const scalar rhoSafe = max(rho, SMALL);
    const scalar TSafe = max(T, SMALL);
    const scalar denom = 1.0/(rhoSafe*rhoSafe)
        + 2.0*bCoeff_/rhoSafe
        - bCoeff_*bCoeff_;
    const scalar alphaT = alpha(TSafe);
    const scalar sqrtAlpha = std::sqrt(alphaT);

    const scalar dPdrhoT =
        R_*TSafe/sqr(1.0 - bCoeff_*rhoSafe)
      - 2.0*aCoeff_*alphaT*rhoSafe*(1.0 + bCoeff_*rhoSafe)
       /sqr(max(denom*rhoSafe*rhoSafe, SMALL));

    const scalar dPdT =
        R_/(1.0/rhoSafe - bCoeff_)
      + aCoeff_*fw_*sqrtAlpha
       /(std::sqrt(max(TSafe*Tc_, SMALL))*max(denom, SMALL));

    return dPdrhoT*gradRho + dPdT*gradT;
}

} // End namespace Foam

// ************************************************************************* //
