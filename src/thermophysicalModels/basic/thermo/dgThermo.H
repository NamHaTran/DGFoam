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
    Foam::dgThermo

Description
    Base manager for thermophysical pipelines in DG:
    - Owns/coordinates: eqnOfState, thermoLaw, transportLaw
    - Provides a common wrapper 'update(...)' taking context inputs/outputs
    - Stores last-evaluated properties: R, Cp, h, e, mu, kappa, Pr, a
      (these can be overridden by derived implementations)

Notes
    - eqnOfState/thermoLaw/transportLaw are your runtime-selectable models.
    - eqnOfState is assumed to expose at least R(). (suy đoán)

SourceFiles
    dgThermo.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgThermo_H
#define Foam_dgThermo_H

#include "word.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "vector.H"
#include "Ostream.H"
#include "dgGeomMesh.H"
#include "eqnOfState.H"
#include "thermoLaw.H"
#include "transportLaw.H"
#include "energy.H"
#include "GaussField.H"

// Forward decl to avoid heavy include
namespace Foam { class dgGeomMesh; }
namespace Foam { class eqnOfState; }
namespace Foam { class thermoLaw; }
namespace Foam { class transportLaw; }

namespace Foam
{

/*------------------------------ Context I/O --------------------------------*/

//- Inputs for dgThermo::update(); pointers are optional.
//  Derived classes check the combination they require.
struct dgThermoInputs
{
    // Conserved variables (per-cell)
    const GaussField<scalar>* rhoC{nullptr};
    const GaussField<vector>* rhoU{nullptr};
    const GaussField<scalar>* rhoE{nullptr};
    // Primitive variables (per-cell)
    const GaussField<scalar>* p{nullptr};
    const GaussField<vector>* U{nullptr};
    const GaussField<scalar>* T{nullptr};
};

//- Outputs from dgThermo::update()
struct dgThermoOutputs
{
    // Primitive (optional write-back by pipelines)
    GaussField<scalar>* rho{nullptr};
    GaussField<vector>* U{nullptr};
    GaussField<scalar>* p{nullptr};
    GaussField<scalar>* T{nullptr};

    // Thermo/transport outputs (required)
    GaussField<scalar>* a{nullptr};
    GaussField<scalar>* mu{nullptr};
    GaussField<scalar>* kappa{nullptr};
    GaussField<scalar>* Cp{nullptr};
    GaussField<scalar>* h{nullptr};
    GaussField<scalar>* e{nullptr};
    GaussField<scalar>* Pr{nullptr};
    GaussField<scalar>* gamma{nullptr};
};

/*---------------------------------------------------------------------------*\
                             Class dgThermo Declaration
\*---------------------------------------------------------------------------*/

class dgThermo
{
protected:

    // Attributes

        //- Model type name (for logging/debug)
        word name_;

        //- Root dictionary (immutable snapshot)
        const dictionary dict_;

        //- Reference to mesh (non-owning)
        dgGeomMesh& mesh_;

        //- Owned models
        autoPtr<eqnOfState> eqnState_;
        autoPtr<thermoLaw>       thermo_;
        autoPtr<transportLaw>    transport_;
        autoPtr<energy>     energy_;

        //- Cached properties (last update)
        scalar R_;                   // Gas constant
        GaussField<scalar> Cp_;      // Specific heat Cp
        GaussField<scalar> h_;       // Enthalpy
        GaussField<scalar> e_;       // Internal energy
        GaussField<scalar> mu_;      // Dynamic viscosity
        GaussField<scalar> kappa_;   // Thermal conductivity
        GaussField<scalar> Pr_;      // Prandtl number
        GaussField<scalar> a_;       // Speed of sound
        GaussField<scalar> gamma_;   // Heat capacity ratio

public:

    //- Runtime type info
    TypeName("dgThermo");

    //- Runtime selection table
    declareRunTimeSelectionTable
    (
        autoPtr, dgThermo, dictionary,
        (const word& name, const dictionary& dict, dgGeomMesh& mesh),
        (name, dict, mesh)
    );

    //- Construct from type name, dict, mesh
    dgThermo(const word& name, const dictionary& dict, dgGeomMesh& mesh);

    //- Disallow default construction
    dgThermo() = delete;

    //- Virtual destructor
    virtual ~dgThermo() = default;

    //- Factory
    static autoPtr<dgThermo> New
    (
        const word& name,
        const dictionary& dict,
        dgGeomMesh& mesh
    );

    // Setup / model wiring (derived will override to construct models)

        //- Initialize inner models (eqnOfState, thermoLaw, transportLaw)
        //  Base does nothing; derived must implement.
        virtual void initModels() = 0;

        //- Optional: validate allowed combinations (derived)
        virtual void validateModelCombo() {}

    // Core API

        //- Wrapper to compute properties from inputs and write into outputs
        //  Derived classes implement the pipeline and must fill outputs.
        virtual void update
        (
            const dgThermoInputs& in,
            dgThermoOutputs& out
        ) = 0;

    // Accessors (can be overridden)

        inline const eqnOfState&      eos() const { return eqnState_(); }
        inline const thermoLaw&       thermo() const { return thermo_(); }
        inline const transportLaw&    transport() const { return transport_(); }

        inline scalar R()                 const { return R_; }
        inline GaussField<scalar> Cp()    const { return Cp_; }
        inline GaussField<scalar> h()     const { return h_; }
        inline GaussField<scalar> e()     const { return e_; }
        inline GaussField<scalar> mu()    const { return mu_; }
        inline GaussField<scalar> kappa() const { return kappa_; }
        inline GaussField<scalar> Pr()    const { return Pr_; }
        inline GaussField<scalar> a()     const { return a_; }
        inline GaussField<scalar> gamma() const { return gamma_; }

        // For derived classes to set cached values conveniently
        inline void setR(const scalar v){ R_=v; }
        inline void setCp(const scalar v){ Cp_=v; }
        inline void setH(const scalar v){ h_=v; }
        inline void setE(const scalar v){ e_=v; }
        inline void setMu(const scalar v){ mu_=v; }
        inline void setKappa(const scalar v){ kappa_=v; }
        inline void setPr(const scalar v){ Pr_=v; }
        inline void setA(const scalar v){ a_=v; }
        inline void setGamma(const scalar v){ gamma_=v; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
