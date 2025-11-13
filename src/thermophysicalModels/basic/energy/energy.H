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
    Foam::energy

Description
    Base class for energy model.

SourceFiles
    energy.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_energy_H
#define Foam_energy_H

#include "scalar.H"
#include "word.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class energy Declaration
\*---------------------------------------------------------------------------*/

class energy
{
protected:

    //- Name of the model (for debugging/logging)
    word name_;

    //- Dictionary for model coefficients (copied at construction)
    const dictionary coeffs_;

public:

    //- Runtime type information
    TypeName("energy");

    //- Declare runtime constructor table
    declareRunTimeSelectionTable
    (
        autoPtr, energy, dictionary,
        (const word& name, const dictionary& dict),
        (name, dict)
    );

    //- Construct from name and dictionary
    energy(const word& name, const dictionary& dict);

    //- Destructor
    virtual ~energy() = default;

    //- Selector: create a new energy model
    static autoPtr<energy> New(const word& name, const dictionary& dict);

    //- Return temperature T from energy (e or h), with Cp and Cv passed in
    virtual GaussField<scalar> T
    (
        const GaussField<scalar>& he,
        const GaussField<scalar>& Cv,
        const GaussField<scalar>& Cp
    ) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //

