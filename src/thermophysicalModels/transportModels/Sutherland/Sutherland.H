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
    Classical Sutherland-law transport model:
    μ(T)  = As * T^(3/2) / (T + S)
    Pr    = constant
    κ(T)  = μ Cp / Pr      (classical Navier–Stokes)

    All evaluations are per-cell GaussField operations.

SourceFiles
    Sutherland.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Sutherland_H
#define Foam_Sutherland_H

#include "transportLaw.H"
#include "word.H"
#include "scalar.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Sutherland Declaration
\*---------------------------------------------------------------------------*/

class Sutherland
:
    public transportLaw
{
private:

    //- Sutherland coefficient [kg/(m·s·sqrt(K))]
    scalar As_;

    //- Sutherland temperature [K]
    scalar S_;

    //- Constant Prandtl number [-]
    scalar Pr0_;

public:

    //- Runtime type information
    TypeName("Sutherland");

    //- Construct from model name and dictionary
    Sutherland
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    //- Destructor
    virtual ~Sutherland() = default;

    //- μ(T)
    virtual void calcMu
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& mu
    ) const override;

    //- κ(T)
    virtual void calcKappa
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& kappa
    ) const override;

    //- Pr(T)
    virtual void calcPr
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& Pr
    ) const override;

    virtual void read() override;

    //- Set group
    virtual bool isSutherlandTransportLaw() const override { return true; }
};

} // End namespace Foam

#endif


// ************************************************************************* //
