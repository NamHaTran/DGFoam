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

Class
    Foam::powerVHS

Description
    Variable-Hard-Sphere (power-law) viscosity model:
      - mu(T) = muRef * (T/TRef)^omega
      - muRef = 15*sqrt(pi*m*kB*TRef) / [ 2*pi*dRef^2*(5-2*omega)*(7-2*omega) ]
      - Pr(T) = Pr0 (constant)
      - kappa(Cp,T) is computed in base transportLaw: mu(T)*Cp/Pr(T)

    Notes:
      - molMass is in g/mol (converted internally).
      - dRef in [m]; TRef in [K]; kB in [J/K]; NA in [1/mol].

    Example:
    transportCoeffs
    {
        type     powerVHS;

        // Required
        molMass  28.0134;     // [g/mol] N2
        dRef     3.7e-10;     // [m]
        TRef     300;         // [K]
        omega    0.74;        // [-]
        Pr       0.72;        // [-]

        // Optional overrides
        // kB    1.380649e-23;   // [J/K]
        // NA    6.02214076e23;  // [1/mol]
    }

SourceFiles
    powerVHS.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_powerVHS_H
#define Foam_powerVHS_H

#include "transportLaw.H"
#include "word.H"
#include "scalar.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class powerVHS Declaration
\*---------------------------------------------------------------------------*/

class powerVHS : public transportLaw
{
private:

    // Coefficients (from dictionary)

        //- Molecular mass [g/mol] (converted to kg/molecule internally)
        scalar molMass_;

        //- Reference hard-sphere diameter [m]
        scalar dRef_;

        //- Reference temperature [K]
        scalar TRef_;

        //- Power-law exponent omega [-]
        scalar omega_;

        //- Constant Prandtl number [-]
        scalar Pr0_;

        //- Boltzmann constant [J/K] (optional override in dict)
        scalar kB_;

        //- Avogadro constant [1/mol] (optional override in dict)
        scalar NA_;

    // Derived/precomputed

        //- Reference dynamic viscosity [kg/(m·s)] at TRef_
        scalar muRef_;

public:

    //- Runtime type information
    TypeName("powerVHS");

    // Constructors

        //- Construct from model name and dictionary; parse coefficients
        powerVHS(const word& name, const dictionary& dict);

    //- Destructor
    ~powerVHS() override = default;

    // Member Functions

        //- Dynamic viscosity mu(T) = muRef * (T/TRef)^omega
        scalar mu(const scalar T) const override;

        //- Prandtl number (constant)
        scalar Pr(const scalar T) const override;

        //- Read/parse coeff() and precompute muRef_
        void read() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
