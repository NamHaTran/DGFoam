/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::dgThermoConservative

Description
    Thermodynamics wrapper for DG solvers based on conservative variables.

    This model reconstructs all primitive thermodynamic quantities from the
    conserved fields (rho, rhoU, E) at each Gauss point of every DG element.
    It provides temperature, pressure, internal energy, enthalpy, heat capacity
    ratio, speed of sound, and transport properties in a fully cell-local and
    Gauss-point-resolved manner.

    dgThermoConservative is designed for compressible DG formulations using
    the conservative form of the Navier–Stokes equations. The model supports
    ideal-gas and thermal-perfect-gas equations of state, together with
    compatible thermo/transport/energy laws (e.g. perfect-gas Cp(T), kinetic
    VHS/VSS, Sutherland viscosity, power-law transport).

    Typical usage:
      - Reconstruct T, p, e, h, Cp, Cv, gamma, a
        from (rho, rhoU, rhoE)
      - Compute transport coefficients mu, kappa, Pr
      - Provide cell-local thermo support for flux solvers, Riemann solvers,
        limiters, viscous operators, and shock capturing schemes.

    This model is intended for DG compressible solvers and does not support
    primitive-variable formulations or incompressible models.

    Below are availabie mapping between classes of dgThermoConservative
     - IdealGas EOS → PerfectGasThermo / AbsoluteEnthalpy / Kinetic
     - ThermalPerfectGas EOS → KineticThermo
     - RealGas EOS → RealGasThermo
     - Incompressible EOS → (not supported by DG conservative)

     - PerfectGasThermo → Sutherland / Polynomial / Constant
     - AbsoluteEnthalpyThermo → Sutherland / Polynomial / Constant
     - RealGasThermo → RealGasTransportLaw
     - KineticThermo → KineticTransportLaw / Constant

     - All thermoLaw → energyInternal / energyEnthalpy
     (But DG conservative uses only energyInternal)


SourceFiles
    dgThermoConservative.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgThermoConservative_H
#define Foam_dgThermoConservative_H

#include "dictionary.H"
#include "dgGeomMesh.H"
#include "dgField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "word.H"
#include "autoPtr.H"
#include "scalar.H"
#include "vector.H"
#include "Ostream.H"
#include "dgGeomMesh.H"
#include "eqnOfState.H"
#include "thermoLaw.H"
#include "transportLaw.H"
#include "energy.H"

namespace Foam
{

class eqnOfState;
class thermoLaw;
class transportLaw;
class energy;

/*---------------------------------------------------------------------------*\
                       Class dgThermoConservative Declaration
\*---------------------------------------------------------------------------*/

class dgThermoConservative
{
protected:

    //=== Attributes =========================================================//

    //- Model type name
    word name_;

    //- Root dictionary (immutable snapshot)
    const dictionary dict_;

    //- Reference to geometric DG mesh (not owning)
    const dgGeomMesh& mesh_;

    //- Owned EOS model
    autoPtr<eqnOfState> eqnState_;

    //- Owned thermo law
    autoPtr<thermoLaw> thermo_;

    //- Owned transport model
    autoPtr<transportLaw> transport_;

    //- Owned energy model
    autoPtr<energy> energy_;

    //- Gas constant
    scalar R_;

    //- Cp and Cv
    dgField<scalar> Cp_;
    dgField<scalar> Cv_;

    //- Enthalpy
    dgField<scalar> h_;

    //- Internal energy
    dgField<scalar> e_;

    //- Dynamic viscosity
    dgField<scalar> mu_;

    //- Thermal conductivity
    dgField<scalar> kappa_;

    //- Prandtl number
    dgField<scalar> Pr_;

    //- Speed of sound
    dgField<scalar> a_;

    //- Heat capacity ratio
    dgField<scalar> gamma_;

    //- Conservative variables
        //- rho
        const dgField<scalar>& rho_;
        
        //- rhoU
        const dgField<vector>& rhoU_;
        
        //- E
        const dgField<scalar>& E_;
        
        //- p (to store result)
        dgField<scalar>& p_;

        //- T (to store result)
        dgField<scalar>& T_;

public:

    //=== Runtime type information ==========================================//

    TypeName("dgThermoConservative");

    // Declare selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dgThermoConservative,
        dictionary,
        (const word& name, const dictionary& dict, const dgGeomMesh& mesh),
        (name, dict, mesh)
    );

    //- Disallow default construction
    dgThermoConservative() = delete;

    //- Main constructor
    dgThermoConservative
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //- Virtual destructor
    virtual ~dgThermoConservative() = default;

    //- Construct from dictionary
    static autoPtr<dgThermoConservative> New
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );


    //=== Model Setup Interface =============================================//

    //- Construct EOS, thermo, transport, energy models (pure virtual)
    virtual void initModels() = 0;

    //- Optional validation for model combinations
    virtual void validateModels() = 0;


    //=== Core API ===========================================================//

    //- Update cell thermodynamics (pure virtual)
    virtual void update(const label& cellID) = 0;


    //=== Accessors ==========================================================//

    //- EOS, thermo, transport
    inline const eqnOfState& eos() const     { return eqnState_(); }
    inline const thermoLaw& thermo() const   { return thermo_(); }
    inline const transportLaw& transport() const { return transport_(); }

    //- Gas constant
    inline scalar R() const { return R_; }

    /*
    // Helper macro to generate field accessors (ref + const-ref)
    #define dgThermoConservative_FIELD_ACCESSOR(fieldName)                               \
        inline dgField<scalar>& fieldName() { return fieldName##_; }         \
        inline const dgField<scalar>& fieldName() const { return fieldName##_; } \
        inline GaussField<scalar>& fieldName(const label cellID)             \
        { return fieldName##_.gaussFields()[cellID]; }                                     \
        inline const GaussField<scalar>& fieldName(const label cellID) const \
        { return fieldName##_.gaussFields()[cellID]; }

    dgThermoConservative_FIELD_ACCESSOR(Cp)
    dgThermoConservative_FIELD_ACCESSOR(Cv)
    dgThermoConservative_FIELD_ACCESSOR(h)
    dgThermoConservative_FIELD_ACCESSOR(e)
    dgThermoConservative_FIELD_ACCESSOR(mu)
    dgThermoConservative_FIELD_ACCESSOR(kappa)
    dgThermoConservative_FIELD_ACCESSOR(Pr)
    dgThermoConservative_FIELD_ACCESSOR(a)
    dgThermoConservative_FIELD_ACCESSOR(gamma)

    #undef dgThermoConservative_FIELD_ACCESSOR
    */

    inline const dgField<scalar>& Cp() const { return Cp_; }
    inline const dgField<scalar>& Cv() const { return Cv_; }
    inline const dgField<scalar>& h() const { return h_; }
    inline const dgField<scalar>& e() const { return e_; }
    inline const dgField<scalar>& mu() const { return mu_; }
    inline const dgField<scalar>& kappa() const { return kappa_; }
    inline const dgField<scalar>& Pr() const { return Pr_; }
    inline const dgField<scalar>& a() const { return a_; }
    inline const dgField<scalar>& gamma() const { return gamma_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //

