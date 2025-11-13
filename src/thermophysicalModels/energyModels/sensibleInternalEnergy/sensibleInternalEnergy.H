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

class sensibleInternalEnergy : public energy
{
protected:

    //- Optional: save coefficients for using in the future
    dictionary coeffDict_;

public:

    //- Runtime type information
    TypeName("sensibleInternalEnergy");

    //- Constructor
    sensibleInternalEnergy(const word& name, const dictionary& dict);

    //- Override method: compute T from internal energy and Cv
    GaussField<scalar> T(const GaussField<scalar>& he, const GaussField<scalar>& Cv, const GaussField<scalar>& Cp) const override;
};

} // namespace Foam

#endif


