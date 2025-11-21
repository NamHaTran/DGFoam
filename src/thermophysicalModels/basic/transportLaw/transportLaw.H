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
    Transport laws define the molecular transport properties used in the
    DG thermophysical system. They provide consistent evaluation of dynamic
    viscosity μ(T,ρ), thermal conductivity κ(T,ρ), Prandtl number, and
    related transport quantities at Gauss points.

    Transport models are grouped into the following categories:

    1. Constant-transport models
    - μ, κ, and Pr are constant.
    - No dependence on temperature or density.
    - Useful for verification cases, laminar flows, and GPU-optimized
    simulations.
    - Compatible with all thermo and EOS models.

    2. Sutherland-type ideal-gas transport
    - Viscosity follows Sutherland's law:
    μ = μ₀ (T/T₀)^(3/2) (T₀ + S) / (T + S).
    - Thermal conductivity is computed via μ Cp / Pr.
    - Valid for ideal gases at moderate temperatures.
    - Compatible with perfect-gas thermo and ideal-gas or
    thermal-perfect-gas EOS.

    3. Polynomial / power-law transport
    - Transport properties vary with temperature using empirical
    correlations, e.g.:
    μ(T) = a₀ + a₁ T + a₂ T² + ...,
    μ(T) = a T^b.
    - Suitable for high-temperature gases, combustion, and fluids
    where viscosity changes strongly with temperature.
    - Works with perfect-gas, constant-CpCv, absolute-enthalpy,
    and selected real-gas thermo models.

    4. Real-gas transport
    - μ and κ depend on both temperature and density:
    μ = μ(T, ρ),   κ = κ(T, ρ).
    - Required for supercritical fluids, hydrocarbons, cryogenic
    fluids, and high-pressure real gases.
    - Must be used together with real-gas EOS models
    (Peng–Robinson, SRK, RK, etc.).

    5. Kinetic / rarefied-gas transport (VHS/VSS/BGK/ES-BGK)
    - Based on kinetic theory of gases.
    - Viscosity and thermal conductivity follow power-law relations:
    μ ∝ T^ω,
    κ = μ Cp / Pr.
    - Includes molecular parameters (collision diameter, molecular
    mass, Boltzmann constant).
    - Required for rarefied-gas flows, high-altitude aerothermodynamics,
    and BGK-type models.
    - Compatible with kinetic thermo and thermal-perfect-gas EOS.

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
#include "dgGeomMesh.H"
#include "dgMath.H"
#include "thermoLaw.H"

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

    //- Reference to DG mesh (non-owning)
    const dgGeomMesh& mesh_;

    //- Reference to thermoLaw
    const thermoLaw& thermo_;

public:

    //- Runtime type information
    TypeName("transportLaw");

    //- Runtime selection table for dictionary-based construction
    declareRunTimeSelectionTable
    (
        autoPtr, transportLaw, dictionary,
        (
            const word& name, 
            const dictionary& dict, 
            const dgGeomMesh& mesh,
            const thermoLaw& thermo
        ),
        (name, dict, mesh, thermo)
    );

    //- Construct from model name and coefficient dictionary
    //  * Store name_ and coeff_ (no validation here).
    transportLaw
    (
        const word& name, 
        const dictionary& dict, 
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    //- Disallow default construction
    transportLaw() = delete;

    //- Virtual destructor
    virtual ~transportLaw() = default;

    //- Selector: create a transportLaw from name and dictionary
    static autoPtr<transportLaw> New
    (
        const word& name, 
        const dictionary& dict, 
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    // Accessors
        //- Compute viscosity μ(T) [kg/(m·s)]
        virtual void calcMu
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& mu
        ) const = 0;

        //- Compute thermal conductivity κ(T) [W/(m·K)]
        virtual void calcKappa
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& kappa
        ) const = 0;

        //- Compute Prandtl number Pr(T)
        virtual void calcPr
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& Pr
        ) const = 0;

    // I/O
        //- Read coefficients (default = do nothing)
        virtual void read();

    // Transport law groups

        //- Constant transport
        virtual bool isConstantTransportLaw() const { return false; }

        //- Sutherland’s law (ideal-gas viscosity)
        virtual bool isSutherlandTransportLaw() const { return false; }

        //- Polynomial/Power-law transport
        virtual bool isPolynomialTransportLaw() const { return false; }

        //- Real-gas transport
        virtual bool isRealGasTransportLaw() const { return false; }

        //- Kinetic / Rarefied transport
        virtual bool isKineticTransportLaw() const { return false; }

};

} // End namespace Foam

#endif

// ************************************************************************* //
