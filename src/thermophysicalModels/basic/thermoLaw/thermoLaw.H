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
    Base class for thermodynamic laws used in DG thermophysical models.

    This class provides the generic interface for:
      - Converting between temperature and energy forms (e, h).
      - Computing Cp(T), Cv(T), gamma(T), and speed of sound.
      - Computing internal energy and enthalpy at Gauss points.

    Thermodynamic models are grouped into the following categories:

    1. Perfect-gas thermo
       - Cp(T) and Cv(T) depend only on temperature.
       - No formation enthalpy is included.
       - h(T) = ∫ Cp(T) dT,  e(T) = ∫ Cv(T) dT.
       - Compatible with ideal-gas EOS (p = rho*R*T).
       - Used in compressible flows and high-temperature gases.

    2. Absolute-enthalpy thermo
       - Enthalpy includes chemical formation energy:
           h(T) = h_f0 + ∫ Cp(T) dT.
       - Required for combustion and reacting flows.
       - Works with ideal-gas or selected real-gas EOS.

    3. Real-gas thermo
       - Cp, Cv, h, and e depend on both temperature and density.
       - Must be used together with real-gas EOS
         (Peng–Robinson, SRK, Redlich–Kwong, etc.).
       - Required for high-pressure and supercritical fluids.

    4. Kinetic / rarefied-gas thermo (VHS/VSS/BGK/ES-BGK)
       - Based on kinetic theory and molecular degrees of freedom.
       - e = (f/2)*R*T,  h = ((f+2)/2)*R*T.
       - Cp and Cv follow kinetic DOF relations.
       - Compatible with thermal-perfect-gas EOS.
       - Used in rarefied gas and hypersonic applications.

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
#include "GaussField.H"
#include "dgGeomMesh.H"
#include "eqnOfState.H"
#include "dgMath.H"

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

    //- Reference to DG mesh (non-owning)
    const dgGeomMesh& mesh_;

    //- Reference to eqnOfState
    const eqnOfState& eos_;

public:

    //- Runtime type information
    TypeName("thermoLaw");

    //- Runtime selection table for dictionary-based construction
    declareRunTimeSelectionTable
    (
        autoPtr, thermoLaw, dictionary,
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh,
            const eqnOfState& eos
        ),
        (name, dict, mesh, eos)
    );

    //- Construct from model name and coefficient dictionary
    thermoLaw
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const eqnOfState& eos
    );

    //- Virtual destructor
    virtual ~thermoLaw() = default;

    //- Selector
    static autoPtr<thermoLaw> New
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const eqnOfState& eos
    );

    //====================================================================//
    //                      General Thermo Interface
    //
    //   All functions:
    //     - receive cellI
    //     - input fields = const GaussField<scalar>&
    //     - output field = GaussField<scalar>&
    //     - return void (write directly into output)
    //
    //   This interface is compatible with all models:
    //     constantCp, perfectGas, JANAF, vibrational, polyatomic, ...
    //====================================================================//

    //- Calculate temperature T from internal energy e(T) or enthalpy h(T)
    virtual void calcT
    (
        const label cellI,
        const GaussField<scalar>& energy,       // input: e or h depending on model
        GaussField<scalar>& T                   // output: temperature
    ) const = 0;

    //- Calculate Cp(T)
    virtual void calcCp
    (
        const label cellI,
        const GaussField<scalar>& T,            // input: temperature
        GaussField<scalar>& Cp                  // output: Cp(T)
    ) const = 0;

    //- Calculate Cv(T)
    virtual void calcCv
    (
        const label cellI,
        const GaussField<scalar>& T,            // input: temperature
        GaussField<scalar>& Cv                  // output: Cv(T)
    ) const = 0;

    //- Calculate gamma(T) = Cp(T)/Cv(T)
    virtual void calcGamma
    (
        const label cellI,
        const GaussField<scalar>& Cp,
        const GaussField<scalar>& Cv,
        GaussField<scalar>& gamma
    ) const = 0;

    //- Calculate internal energy e(T)
    virtual void calcInternalE
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& e
    ) const = 0;

    //- Calculate enthalpy h(T)
    virtual void calcH
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& h
    ) const = 0;

    //- Calculate speedOfSound
    virtual void calcSpeedOfSound
    (
        const label cellI,
        const GaussField<scalar>& T,
        const GaussField<scalar>& gamma,
        GaussField<scalar>& a
    ) const = 0;

    //- ThermoLaw groups
        // Perfect-gas type: Cp(T), Cv(T) but h=∫Cp dT (no formation energy)
        virtual bool isPerfectGasThermo() const { return false; }

        // Absolute enthalpy model: h includes formation enthalpy
        virtual bool isAbsoluteEnthalpyThermo() const { return false; }

        // Real-gas / high-temperature thermo
        virtual bool isRealGasThermo() const { return false; }

        // Kinetic/rarefied thermo (VHS/VSS)
        virtual bool isKineticThermo() const { return false; }
};

} // End namespace Foam

#endif

// ************************************************************************* //

