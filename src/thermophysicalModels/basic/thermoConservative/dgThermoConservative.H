/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgThermoConservative

Description
    Conservative-variable thermodynamics model for DG solvers.

SourceFiles
    dgThermoConservative.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgThermoConservative_H
#define Foam_dgThermoConservative_H

#include "dictionary.H"
#include "regIOobject.H"
#include "dgGeomMesh.H"
#include "dgField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "word.H"
#include "scalar.H"
#include "vector.H"
#include "Ostream.H"
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

/**
 * \class Foam::dgThermoConservative
 * \brief Conservative thermodynamics wrapper for DG solvers
 *
 * \details
 * dgThermoConservative reconstructs primitive thermodynamic quantities
 * at DG Gauss points from conservative variables:
 *
 * \f[
 *   (\rho,\; \rho\mathbf{U},\; E)
 * \;\longrightarrow\;
 * (T,\; p,\; e,\; h,\; C_p,\; C_v,\; \gamma,\; a)
 * \f]
 *
 * The model operates in a fully cell-local and Gauss-point-resolved
 * manner and provides thermodynamic and transport properties required
 * by compressible DG formulations.
 *
 * Supported model structure:
 * - Equation of state (eqnOfState)
 * - Thermodynamic law (thermoLaw)
 * - Transport law (transportLaw)
 * - Energy formulation (energyInternal only)
 *
 * The class is intended exclusively for conservative-variable DG
 * formulations of compressible flow and does not support primitive
 * or incompressible models.
 */
class dgThermoConservative
:
    public regIOobject
{
protected:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Model instance name */
        word name_;

        /** Root dictionary (immutable snapshot) */
        const dictionary dict_;

        /** Reference to DG geometric mesh (non-owning) */
        const dgGeomMesh& mesh_;

        /** Equation of state model (owning) */
        autoPtr<eqnOfState> eqnState_;

        /** Thermodynamic law model (owning) */
        autoPtr<thermoLaw> thermo_;

        /** Transport law model (owning) */
        autoPtr<transportLaw> transport_;

        /** Energy formulation model (owning) */
        autoPtr<energy> energy_;

        /** Gas constant */
        scalar R_;

        /** Heat capacity at constant pressure */
        dgField<scalar> Cp_;

        /** Heat capacity at constant volume */
        dgField<scalar> Cv_;

        /** Specific enthalpy */
        dgField<scalar> h_;

        /** Specific internal energy */
        dgField<scalar> e_;

        /** Dynamic viscosity */
        dgField<scalar> mu_;

        /** Thermal conductivity */
        dgField<scalar> kappa_;

        /** Prandtl number */
        dgField<scalar> Pr_;

        /** Speed of sound */
        dgField<scalar> a_;

        /** Heat-capacity ratio */
        dgField<scalar> gamma_;

        /** Conservative density field (read-only) */
        const dgField<scalar>& rho_;

        /** Conservative momentum field (read-only) */
        const dgField<vector>& rhoU_;

        /** Conservative total energy field (read-only) */
        const dgField<scalar>& E_;

        /** Pressure field (output storage) */
        dgField<scalar>& p_;

        /** Temperature field (output storage) */
        dgField<scalar>& T_;

public:

    //====================================================================//
    //                     Runtime type information
    //====================================================================//

    TypeName("dgThermoConservative");

    declareRunTimeSelectionTable
    (
        autoPtr,
        dgThermoConservative,
        dictionary,
        (const word& name, const dictionary& dict, const dgGeomMesh& mesh),
        (name, dict, mesh)
    );

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** Disallow default construction */
        dgThermoConservative() = delete;

        /**
         * \brief Construct from name, dictionary, and DG mesh
         *
         * \param name Model instance name
         * \param dict Thermodynamics dictionary
         * \param mesh DG geometric mesh
         */
        dgThermoConservative
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

        /** Virtual destructor */
        virtual ~dgThermoConservative() = default;

        /**
         * \brief Construct via runtime selection
         *
         * \param name Model instance name
         * \param dict Thermodynamics dictionary
         * \param mesh DG geometric mesh
         */
        static autoPtr<dgThermoConservative> New
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

    //====================================================================//
    //                         Model initialization
    //====================================================================//

        /**
         * \brief Construct EOS, thermo, transport, and energy models
         */
        virtual void initModels() = 0;

        /**
         * \brief Validate selected model combinations
         */
        virtual void validateModels() = 0;

    //====================================================================//
    //                              Core API
    //====================================================================//

        /**
         * \brief Update thermodynamic fields for a cell
         *
         * \param cellID DG cell index
         */
        virtual void update(const label& cellID) = 0;

        /**
         * \brief Write model data
         *
         * \param os Output stream
         * \return True if write succeeded
         */
        virtual bool writeData(Ostream& os) const override;

    //====================================================================//
    //                              Accessors
    //====================================================================//

        /** Return equation-of-state model */
        inline const eqnOfState& eos() const { return eqnState_(); }

        /** Return thermodynamic law */
        inline const thermoLaw& thermo() const { return thermo_(); }

        /** Return transport law */
        inline const transportLaw& transport() const { return transport_(); }

        /** Return gas constant */
        inline scalar R() const { return R_; }

        /** Return Cp field */
        inline const dgField<scalar>& Cp() const { return Cp_; }

        /** Return Cv field */
        inline const dgField<scalar>& Cv() const { return Cv_; }

        /** Return enthalpy field */
        inline const dgField<scalar>& h() const { return h_; }

        /** Return internal energy field */
        inline const dgField<scalar>& e() const { return e_; }

        /** Return viscosity field */
        inline const dgField<scalar>& mu() const { return mu_; }

        /** Return thermal conductivity field */
        inline const dgField<scalar>& kappa() const { return kappa_; }

        /** Return Prandtl number field */
        inline const dgField<scalar>& Pr() const { return Pr_; }

        /** Return speed-of-sound field */
        inline const dgField<scalar>& a() const { return a_; }

        /** Return heat-capacity ratio field */
        inline const dgField<scalar>& gamma() const { return gamma_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
