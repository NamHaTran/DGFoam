/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "powerVHS.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "IOstreams.H"
#include <cmath>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(powerVHS, 0);

// Register this model into transportLaw dictionary table
addToRunTimeSelectionTable(transportLaw, powerVHS, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerVHS::powerVHS
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
:
    transportLaw(name, dict, mesh, thermo),
    molMass_(28.0134),   // air-like (N2)
    dRef_(3.7e-10),      // 3.7 Å
    TRef_(300.0),
    omega_(0.74),
    Pr0_(0.72),
    kB_(1.380649e-23),
    NA_(6.02214076e23),
    muRef_(1.8e-5)       // temporary seed (updated in read)
{
    powerVHS::read();
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::powerVHS::calcMu
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& mu
) const
{
    // Guard: T must be positive everywhere in the cell
    if (T <= scalar(0))
    {
        FatalErrorInFunction
            << "Temperature must be strictly positive in powerVHS::calcMu()." << nl
            << exit(FatalError);
    }

    // μ = μRef * (T/TRef)^ω
    GaussField<scalar> theta(cellI, &mesh_);
    theta = T / TRef_;

    GaussField<scalar> thetaW(cellI, &mesh_);
    thetaW = pow(theta, omega_);

    mu = muRef_ * thetaW;
}


void Foam::powerVHS::calcKappa
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& kappa
) const
{
    // classical energy transport: κ = μ Cp / Pr
  
    GaussField<scalar> mu(cellI, &mesh_);
    calcMu(cellI, T, mu);

    GaussField<scalar> Cp(cellI, &mesh_);
    thermo_.calcCp(cellI, T, Cp);

    kappa = (mu * Cp) / Pr0_;
}


void Foam::powerVHS::calcPr
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& Pr
) const
{
    Pr = Pr0_;   // broadcast to all Gauss points
}


// * * * * * * * * * * * * * * * * Read Coeffs  * * * * * * * * * * * * * * //

void Foam::powerVHS::read()
{
    // Parse input coefficients
    if (coeff_.found("molMass")) molMass_ = readScalar(coeff_.lookup("molMass"));
    if (coeff_.found("dRef"))    dRef_    = readScalar(coeff_.lookup("dRef"));
    if (coeff_.found("TRef"))    TRef_    = readScalar(coeff_.lookup("TRef"));
    if (coeff_.found("omega"))   omega_   = readScalar(coeff_.lookup("omega"));
    if (coeff_.found("Pr"))      Pr0_     = readScalar(coeff_.lookup("Pr"));

    if (coeff_.found("kB")) kB_ = readScalar(coeff_.lookup("kB"));
    if (coeff_.found("NA")) NA_ = readScalar(coeff_.lookup("NA"));


    // Validate
    if (molMass_ <= 0)
        FatalErrorInFunction << "molMass must be positive." << exit(FatalError);

    if (dRef_ <= 0)
        FatalErrorInFunction << "dRef must be positive." << exit(FatalError);

    if (TRef_ <= 0)
        FatalErrorInFunction << "TRef must be positive." << exit(FatalError);

    if (Pr0_ <= 0)
        FatalErrorInFunction << "Pr must be positive." << exit(FatalError);

    // (5 - 2ω)(7 - 2ω) must not be zero
    const scalar a = 5 - 2*omega_;
    const scalar b = 7 - 2*omega_;
    const scalar denomAB = a*b;

    if (mag(denomAB) <= SMALL)
    {
        FatalErrorInFunction
            << "(5 - 2*omega)*(7 - 2*omega) is near zero." << nl
            << exit(FatalError);
    }

    // Compute molecular mass in kg/molecule
    const scalar m = (molMass_ * 1e-3) / NA_;

    // Compute muRef via VHS formulation
    const scalar pi = constant::mathematical::pi;

    const scalar num =
        15.0 * std::sqrt(pi * m * kB_ * TRef_);

    const scalar denom =
        2.0 * pi * dRef_ * dRef_ * denomAB;

    if (mag(denom) <= SMALL)
    {
        FatalErrorInFunction
            << "Denominator in muRef formula is near zero." << nl
            << exit(FatalError);
    }

    muRef_ = num / denom;

    if (muRef_ <= 0)
    {
        FatalErrorInFunction
            << "Computed muRef is non-positive: " << muRef_ << nl
            << exit(FatalError);
    }
}

// ************************************************************************* //

} // End namespace Foam
