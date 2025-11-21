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

    //- Gas constant R [J/(kg.K)]
    scalar R_;

public:

    //- Runtime type information
    TypeName("idealGas");

    //- Runtime selection table
    declareRunTimeSelectionTable
    (
        autoPtr, idealGas, dictionary,
        (const word& name, const dictionary& dict, const dgGeomMesh& mesh),
        (name, dict, mesh)
    );

    //- Construct from dictionary
    idealGas
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //- Destructor
    virtual ~idealGas() = default;

    //- Return gas constant R
    inline scalar R() const override
    {
        return R_;
    }

    //==============================================================//
    //         Implement ideal-gas specific thermodynamic relations
    //==============================================================//

    //- rho = p / (R * T)
    virtual void calcRhoFromPT
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& T,
        GaussField<scalar>& rho
    ) const override;

    //- p = rho * R * T
    virtual void calcPFromRhoT
    (
        const label cellI,
        const GaussField<scalar>& rho,
        const GaussField<scalar>& T,
        GaussField<scalar>& p
    ) const override;

    //- T = p / (rho * R)
    virtual void calcTFromPRho
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& rho,
        GaussField<scalar>& T
    ) const override;

    //==============================================================//
    //              Capability flags (ideal gas supports all)
    //==============================================================//

    virtual bool canCalcRhoFromPT()  const override { return true; }
    virtual bool canCalcPFromRhoT()  const override { return true; }
    virtual bool canCalcTFromPRho()  const override { return true; }

    bool isIdealGas() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
