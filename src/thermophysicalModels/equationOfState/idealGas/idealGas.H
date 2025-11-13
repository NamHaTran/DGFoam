/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::idealGas

Description
    Ideal gas equation of state:
        p = rho * R * T

SourceFiles
    idealGas.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_idealGas_H
#define Foam_idealGas_H

#include "eqnOfState.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class idealGas Declaration
\*---------------------------------------------------------------------------*/

class idealGas
:
    public eqnOfState
{
private:

    //- Specific gas constant [J/(kg.K)]
    scalar R_;

public:

    //- Runtime type information
    TypeName("idealGas");

    //- Construct from name and dictionary
    idealGas(const word& name, const dictionary& dict);

    //- Destructor
    virtual ~idealGas() = default;

    // Member Functions

        //- Return density [kg/m3] from T [K] and p [Pa]
        virtual GaussField<scalar> rho(GaussField<scalar>& T, GaussField<scalar>& p) const override;

        //- Return pressure [Pa] from rho [kg/m3] and T [K]
        virtual GaussField<scalar> p(GaussField<scalar>& rho, GaussField<scalar>& T) const override;

        //- Return temperature [K] from rho [kg/m3] and p [Pa]
        virtual GaussField<scalar> T(GaussField<scalar>& rho, GaussField<scalar>& p) const override;

        //- Return gas constant [J/(kg.K)]
        virtual scalar R() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
