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
    Foam::eqnOfState

Description
    Base class for Equation of State (single-phase).
    - Constructible from a coefficients dictionary.
    - Provides rho(T,p) and p(T,rho); a(T) optional.

SourceFiles
    eqnOfState.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_eqnOfState_H
#define Foam_eqnOfState_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class eqnOfState Declaration
\*---------------------------------------------------------------------------*/

class eqnOfState
{
protected:

    //- Model name (for logging/debug)
    word name_;

    //- Coefficients dictionary (owned/const view)
    const dictionary dict_;

public:

    //- Runtime type information
    TypeName("eqnOfState");

    //- Runtime selection: construct from dictionary
    declareRunTimeSelectionTable
    (
        autoPtr, eqnOfState, dictionary,
        (const word& name, const dictionary& dict),
        (name, dict)
    );

    //- Construct from model name and coefficients dictionary
    eqnOfState(const word& name, const dictionary& dict);

    //- Virtual destructor
    virtual ~eqnOfState() = default;

    //- Factory: create from selection table
    static autoPtr<eqnOfState> New(const word& name, const dictionary& dict);

    //- Density rho(T,p)
    virtual GaussField<scalar> rho(const GaussField<scalar>& T, const GaussField<scalar>& p) const = 0;

    //- Pressure p(T,rho)
    virtual GaussField<scalar> p(const GaussField<scalar>& rho, const GaussField<scalar>& T) const = 0;

    //- Temperature T(p,rho)
    virtual GaussField<scalar> T(const GaussField<scalar>& rho, const GaussField<scalar>& p) const = 0;

    //- Gas cpnstant R
    virtual scalar R() const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
