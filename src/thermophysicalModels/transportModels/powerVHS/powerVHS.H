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
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class powerVHS Declaration
\*---------------------------------------------------------------------------*/

class powerVHS
:
    public transportLaw
{
private:

    // Coefficients read from dictionary

        //- Molecular mass [g/mol] (converted internally to kg/molecule)
        scalar molMass_;

        //- Reference diameter [m]
        scalar dRef_;

        //- Reference temperature [K]
        scalar TRef_;

        //- Power-law exponent ω [-]
        scalar omega_;

        //- Constant Prandtl number [-]
        scalar Pr0_;

        //- Boltzmann constant [J/K] (optional override)
        scalar kB_;

        //- Avogadro constant [1/mol] (optional override)
        scalar NA_;

    // Derived values

        //- Reference viscosity at TRef [kg/(m·s)]
        scalar muRef_;

public:

    //- Runtime type information
    TypeName("powerVHS");

    //- Construct from dictionary (reads + precomputes)
    powerVHS
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    //- Destructor
    ~powerVHS() override = default;


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

    //- Dynamic viscosity μ(T)
    virtual void calcMu
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& mu
    ) const override;

    //- Thermal conductivity κ(T)
    virtual void calcKappa
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& kappa
    ) const override;

    //- Prandtl number (constant)
    virtual void calcPr
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& Pr
    ) const override;

    //- Read coefficients + compute μRef_
    virtual void read() override;

    //- Set group
    virtual bool isKineticTransportLaw() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
