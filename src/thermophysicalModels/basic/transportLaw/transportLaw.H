/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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
    Foam::transportLaw

Description
    Abstract base class for laminar transport properties.
    - Provides interface for dynamic viscosity mu(T) and thermal conductivity kappa(T).
    - Supplies a convenience Prandtl-number evaluation given Cp and T.

    Notes (assumptions - may be adapted for your dict layout):
    - Coefficients are stored under a sub-dictionary "transport" of the provided dict.
    - Derived classes read their own coefficients from coeffs_.

SourceFiles
    transportLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_transportLaw_H
#define Foam_transportLaw_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class transportLaw Declaration
\*---------------------------------------------------------------------------*/

class transportLaw
{
protected:

    //- Name of the model (for logging/debugging)
    word name_;

    //- Coefficient dictionary snapshot (immutable after construction)
    //  Always set from input dict; never empty by design.
    const dictionary coeff_;

public:

    //- Runtime type information
    TypeName("transportLaw");

    //- Runtime selection table for dictionary-based construction
    declareRunTimeSelectionTable
    (
        autoPtr, transportLaw, dictionary,
        (const word& name, const dictionary& dict),
        (name, dict)
    );

    //- Construct from model name and coefficient dictionary
    //  * Store name_ and coeff_ (no validation here).
    transportLaw(const word& name, const dictionary& dict);

    //- Disallow default construction
    transportLaw() = delete;

    //- Virtual destructor
    virtual ~transportLaw() = default;

    //- Selector: create a transportLaw from name and dictionary
    static autoPtr<transportLaw> New(const word& name, const dictionary& dict);

    // Accessors
    //- Dynamic viscosity mu(T) [kg/(m·s)]
    virtual GaussField<scalar> mu(const GaussField<scalar>& T) const = 0;

    //- Thermal conductivity kappa(T) [W/(m·K)]
    virtual GaussField<scalar> kappa(const GaussField<scalar>& Cp, const GaussField<scalar>& T) const;

    //- Prandtl number
    virtual GaussField<scalar> Pr(const GaussField<scalar>& T) const = 0;

    // Reading / I/O

    //- Wrapper for derived classes to parse coeff() into internal parameters
    //  * Base implementation does nothing.
    virtual void read();
};

} // End namespace Foam

#endif

// ************************************************************************* //
