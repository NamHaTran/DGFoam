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
    Foam::thermoLaw

Description
    Base class cho các mô hình nhiệt động (single-phase).
    - Khởi tạo từ dictionary hệ số.
    - Cung cấp giao diện Cp(T), h(T), e(T).
    - Dùng runtime selection table để chọn model lúc chạy.

SourceFiles
    thermoLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_thermoLaw_H
#define Foam_thermoLaw_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class thermoLaw Declaration
\*---------------------------------------------------------------------------*/

class thermoLaw
{
protected:

    //- Name of the model (for logging/debugging)
    word name_;

    //- Coefficient dictionary (read-only after construction)
    const dictionary dict_;

public:

    //- Runtime type information
    TypeName("thermoLaw");

    //- Runtime selection table for dictionary-based construction
    declareRunTimeSelectionTable
    (
        autoPtr, thermoLaw, dictionary,
        (const word& name, const dictionary& dict),
        (name, dict)
    );

    //- Construct from model name and coefficient dictionary
    thermoLaw(const word& name, const dictionary& dict);

    //- Virtual destructor
    virtual ~thermoLaw() = default;

    //- Selector: create a thermoLaw from name and dictionary
    static autoPtr<thermoLaw> New(const word& name, const dictionary& dict);

    //- Specific heat capacity Cp(T) [J/(kg·K)]
    virtual scalar Cp(const scalar T) const = 0;

    //- Sensible enthalpy h(T) [J/kg]
    virtual scalar h(const scalar T) const = 0;

    //- Sensible internal energy e(T) [J/kg]
    virtual scalar e(const scalar T) const = 0;

    //- T from sensible internal energy
    virtual scalar TfrE(const scalar e) const = 0;

    //- T from sensible enthalpy energy
    virtual scalar TfrH(const scalar e) const = 0;

    //- Speed of sound a(T) (optional; default aborts)
    virtual scalar a(const scalar T) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
