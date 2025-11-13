/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright
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

#include "Sutherland.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H"
#include <cmath>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Sutherland, 0);

// Register derived into base transportLaw dictionary-ctor table
addToRunTimeSelectionTable(transportLaw, Sutherland, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from name/dict: set safe defaults then parse coeff_
Foam::Sutherland::Sutherland
(
    const word& name,
    const dictionary& dict
)
:
    transportLaw(name, dict),
    As_(1.458e-5),      // default (air-like) [kg/(m·s)]
    S_(110.4),          // default [K]
    Pr0_(0.72)          // default [-]
{
    // Parse and validate dictionary coefficients
    Foam::Sutherland::read();
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

// mu(T) = muRef * (T/TRef)^(3/2) * (TRef + S)/(T + S)
Foam::GaussField<scalar> Foam::Sutherland::mu(const GaussField<scalar>& T) const
{
    // Guard: T must be positive
    if (T <= scalar(0))
    {
        FatalErrorInFunction
            << "Temperature must be positive. T=" << T << nl
            << exit(FatalError);
    }

    const GaussField<scalar> denom = T + S_;

    // Guard: avoid division by zero
    if (denom == scalar(0))
    {
        FatalErrorInFunction
            << "Denominator (T + S) equals zero. T=" << T << ", S=" << S_ << nl
            << exit(FatalError);
    }

    const GaussField<scalar> muT = As_ * pow(T, scalar(1.5)) / denom;

    return muT;
}


// Pr(T): constant Prandtl number
Foam::GaussField<scalar> Foam::Sutherland::Pr(const GaussField<scalar>& T) const
{
    return Pr0_;
}


// Read coefficients from coeff_ and validate
void Foam::Sutherland::read()
{
    // Expected keys in coeff_:
    if (coeff_.found("As"))    { As_    = readScalar(coeff_.lookup("As"));    }
    if (coeff_.found("S"))     { S_     = readScalar(coeff_.lookup("S"));     }
    if (coeff_.found("Pr"))    { Pr0_   = readScalar(coeff_.lookup("Pr"));    }

    // Basic validation
    if (As_ <= scalar(0))
    {
        FatalErrorInFunction
            << "As must be positive. As=" << As_ << nl
            << exit(FatalError);
    }
    if (S_ <= scalar(0))
    {
        FatalErrorInFunction
            << "S must be positive. S=" << S_ << nl
            << exit(FatalError);
    }
    if (Pr0_ <= scalar(0))
    {
        FatalErrorInFunction
            << "Pr must be positive. Pr=" << Pr0_ << nl
            << exit(FatalError);
    }
}

// ************************************************************************* //

} // End namespace Foam
