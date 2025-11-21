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
    Foam::sensibleInternalEnergy

Description
    sensibleInternalEnergy energy model.
    Energy field is sensible internal energy e = Cv*T

SourceFiles
    sensibleInternalEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sensibleInternalEnergy_H
#define Foam_sensibleInternalEnergy_H

#include "energy.H"
#include "word.H"
#include "dictionary.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sensibleInternalEnergy Declaration
\*---------------------------------------------------------------------------*/

class sensibleInternalEnergy
:
    public energy
{
protected:

    //- Optional coefficient storage
    dictionary coeffDict_;

public:

    //- Runtime type information
    TypeName("sensibleInternalEnergy");

    //- Construct from dictionary
    sensibleInternalEnergy
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    // * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    //- Compute enthalpy h(T)
    virtual void calcEnthalpy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& h
    ) const override;

    //- Compute internal energy e(T)
    virtual void calcEnergy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& e
    ) const override;

    //- Compute T(e)
    virtual void calcTfromEnergy
    (
        const label cellI,
        const GaussField<scalar>& e,
        GaussField<scalar>& T
    ) const override;

    //- Set energy group
    virtual bool energyInternal() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //

