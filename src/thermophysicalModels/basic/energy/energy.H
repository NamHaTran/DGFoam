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
    Energy models in DGFoam determine the primary thermodynamic variable
    evolved by the solver (internal energy or enthalpy) and provide
    consistent conversions between temperature and energy forms.

    The energy system is divided into the following categories:

    1. Internal-energy based models
    - The conserved variable is the sensible internal energy e(T).
    - Temperature is obtained from e(T) using thermoLaw relations.
    - Enthalpy is computed as:  h = e + p/rho.
    - No formation enthalpy is included.
    - Suitable for compressible flows, shock-capturing, DG solvers,
    and real-gas thermodynamics.
    - Compatible with perfect-gas, constant-CpCv, real-gas, and kinetic
    thermo models.

    2. Enthalpy-based models
    - The conserved variable is the sensible enthalpy h(T).
    - Temperature is obtained from h(T).
    - Internal energy is computed as:  e = h - p/rho.
    - Required when the transport or chemistry formulation is expressed
    naturally in terms of enthalpy.
    - Compatible with all thermo categories, including
    absolute-enthalpy thermo.

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
#include "dgGeomMesh.H"
#include "dgMath.H"
#include "thermoLaw.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class energy Declaration
\*---------------------------------------------------------------------------*/

class energy
{
protected:

    //- Model name (for logging/debug)
    word name_;

    //- Coefficient dictionary snapshot
    const dictionary coeffs_;

    //- Reference to DG geometric mesh (non-owning)
    const dgGeomMesh& mesh_;

    //- Reference to thermoLaw (non-owning)
    const thermoLaw& thermo_;

public:

    //- Runtime type information
    TypeName("energy");

    //- Runtime selection table
    declareRunTimeSelectionTable
    (
        autoPtr, energy, dictionary,
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh,
            const thermoLaw& thermo
        ),
        (name, dict, mesh, thermo)
    );

    //- Construct from name, dictionary, mesh, and thermoLaw
    energy
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    //- Destructor
    virtual ~energy() = default;

    //- Selector
    static autoPtr<energy> New
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    // * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

    //- Compute enthalpy from temperature (or vice versa)
    //  Output h(cellI) is written directly
    virtual void calcEnthalpy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& h
    ) const = 0;

    //- Compute energy (e or h depending on model)
    //  Output E(cellI) is written directly
    virtual void calcEnergy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& E
    ) const = 0;

    //- Compute temperature from energy (inversion)
    virtual void calcTfromEnergy
    (
        const label cellI,
        const GaussField<scalar>& E,
        GaussField<scalar>& T
    ) const = 0;

    // Energy groups

        //- Internal energy (based on e)
        virtual bool energyInternal() const { return false; }

        //- Enthalpy energy (based on h)
        virtual bool energyEnthalpy() const { return false; }
};

} // End namespace Foam

#endif

// ************************************************************************* //

