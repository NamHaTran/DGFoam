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

// Register this derived model into transportLaw's dictionary-ctor table
addToRunTimeSelectionTable(transportLaw, powerVHS, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from name/dict, set safe defaults, then parse and precompute
Foam::powerVHS::powerVHS
(
    const word& name,
    const dictionary& dict
)
:
    transportLaw(name, dict),
    // Default coefficients (may be overridden in dict)
    molMass_(28.0134),      // N2 ~ 28.0134 g/mol
    dRef_(3.7e-10),         // m (ví dụ 3.7 Å)
    TRef_(300.0),           // K
    omega_(0.74),           // typically 0.74 for air-like gases
    Pr0_(0.72),             // constant Prandtl number
    kB_(1.380649e-23),      // J/K
    NA_(6.02214076e23),     // 1/mol
    muRef_(1.8e-5)          // seed; will be computed in read()
{
    Foam::powerVHS::read(); // parse coeffs + compute muRef_
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

// Compute mu(T) using power law
Foam::GaussField<scalar> Foam::powerVHS::mu(const GaussField<scalar>& T) const
{
    // Step: basic guard for temperature domain
    if (T <= scalar(0))
    {
        FatalErrorInFunction
            << "Temperature must be positive. T=" << T << nl
            << exit(FatalError);
    }

    // Step: compute mu(T) = muRef_ * (T/TRef_)^omega_
    const GaussField<scalar> theta = T/TRef_;
    return muRef_ * pow(theta, omega_);
}


// Return constant Prandtl number
Foam::GaussField<scalar> Foam::powerVHS::Pr(const GaussField<scalar>& T) const
{
    GaussField<scalar> Pr0(T.cellID(), T.dgMesh(), Pr0_);
    return Pr0;
}


// Parse coefficients and precompute muRef_
void Foam::powerVHS::read()
{
    // * * * Parse user coefficients from coeff_ * * *
    if (coeff_.found("molMass")) { molMass_ = readScalar(coeff_.lookup("molMass")); } // [g/mol]
    if (coeff_.found("dRef"))    { dRef_    = readScalar(coeff_.lookup("dRef"));    } // [m]
    if (coeff_.found("TRef"))    { TRef_    = readScalar(coeff_.lookup("TRef"));    } // [K]
    if (coeff_.found("omega"))   { omega_   = readScalar(coeff_.lookup("omega"));   } // [-]
    if (coeff_.found("Pr"))      { Pr0_     = readScalar(coeff_.lookup("Pr"));      } // [-]

    // Optional overrides of constants
    if (coeff_.found("kB")) { kB_ = readScalar(coeff_.lookup("kB")); } // [J/K]
    if (coeff_.found("NA")) { NA_ = readScalar(coeff_.lookup("NA")); } // [1/mol]

    // * * * Validate inputs * * *
    if (molMass_ <= scalar(0))
    {
        FatalErrorInFunction << "molMass must be positive [g/mol]. molMass=" << molMass_ << nl
                             << exit(FatalError);
    }
    if (dRef_ <= scalar(0))
    {
        FatalErrorInFunction << "dRef must be positive [m]. dRef=" << dRef_ << nl
                             << exit(FatalError);
    }
    if (TRef_ <= scalar(0))
    {
        FatalErrorInFunction << "TRef must be positive [K]. TRef=" << TRef_ << nl
                             << exit(FatalError);
    }
    if (Pr0_ <= scalar(0))
    {
        FatalErrorInFunction << "Pr must be positive. Pr=" << Pr0_ << nl
                             << exit(FatalError);
    }

    // Denominator safety: (5 - 2*omega)*(7 - 2*omega) must be non-zero
    const scalar a = (5 - 2*omega_);
    const scalar b = (7 - 2*omega_);
    const scalar denomAB = a*b;
    if (mag(denomAB) <= SMALL)
    {
        FatalErrorInFunction
            << "(5 - 2*omega)*(7 - 2*omega) is near zero. omega=" << omega_ << nl
            << exit(FatalError);
    }

    // * * * Precompute muRef_ * * *
    // Step 1: convert molMass_ [g/mol] -> m [kg/molecule]
    //         m = (molMass_*1e-3 [kg/mol]) / NA_ [1/mol]
    const scalar m = (molMass_ * 1e-3) / NA_;

    // Step 2: compute muRef with VHS formula
    const scalar pi = constant::mathematical::pi;
    const scalar num = 15.0 * std::sqrt(pi * m * kB_ * TRef_);
    const scalar denom = 2.0 * pi * dRef_ * dRef_ * denomAB;

    if (mag(denom) <= SMALL)
    {
        FatalErrorInFunction
            << "Denominator in muRef formula is near zero. dRef=" << dRef_
            << ", omega=" << omega_ << nl
            << exit(FatalError);
    }

    muRef_ = num/denom;

    // Optional: sanity check
    if (muRef_ <= scalar(0))
    {
        FatalErrorInFunction
            << "Computed muRef is non-positive. muRef=" << muRef_ << nl
            << exit(FatalError);
    }
}

// ************************************************************************* //

} // End namespace Foam
