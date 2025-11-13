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
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Sutherland

Description
    Transport model using Sutherland's law:
      mu(T) = muRef * (T/TRef)^(3/2) * (TRef + S)/(T + S)

    - Pr(T) is treated as constant: Pr0.
    - kappa(Cp,T) is computed by the base class:
        kappa = mu(T) * Cp / Pr(T).

SourceFiles
    Sutherland.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Sutherland_H
#define Foam_Sutherland_H

#include "transportLaw.H"
#include "word.H"
#include "scalar.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Sutherland Declaration
\*---------------------------------------------------------------------------*/

class Sutherland : public transportLaw
{
private:

    // Attributes

        //- Sutherland constant [kg/(m*s*sqrt(K))]
        scalar As_;

        //- Sutherland constant [K]
        scalar S_;

        //- Constant Prandtl number [-]
        scalar Pr0_;

public:

    //- Runtime type information
    TypeName("Sutherland");

    // Constructors

        //- Construct from model name and dictionary; parse coefficients
        Sutherland(const word& name, const dictionary& dict);

    //- Destructor
    ~Sutherland() override = default;

    // Member Functions

        //- Dynamic viscosity mu(T) via Sutherland correlation
        GaussField<scalar> mu(const GaussField<scalar>& T) const override;

        //- Prandtl number (constant)
        GaussField<scalar> Pr(const GaussField<scalar>& T) const override;

        //- Read/parse coefficients from coeff() and validate
        void read() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
