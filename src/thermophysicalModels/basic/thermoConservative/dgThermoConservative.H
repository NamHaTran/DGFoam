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
 * (T,\; p,\; he,\; a)
 * \f]
 *
 * The model operates in a fully cell-local and Gauss-point-resolved
 * manner and provides the stored thermodynamic and transport quantities
 * required by compressible DG formulations. Additional properties such as
 * \f$C_p\f$, \f$C_v\f$, \f$\gamma\f$, \f$\kappa\f$, \f$\mu\f$, and
 * \f$\mathrm{Pr}\f$ are intended to be evaluated on demand through the
 * thermoLaw and transportLaw APIs instead of being stored as dgField data.
 * The stored thermodynamic energy variable \c he may represent either
 * sensible internal energy or sensible enthalpy, depending on the selected
 * energy model.
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

        /** Stored thermodynamic energy variable */
        dgField<scalar> he_;

        /** Speed of sound */
        dgField<scalar> a_;

        /** Conservative density field (read-only) */
        const dgField<scalar>& rho_;

        /** Conservative momentum field (read-only) */
        const dgField<vector>& rhoU_;

        /** Conservative total energy field (read-only) */
        const dgField<scalar>& E_;

        /** Velocity field (output storage) */
        dgField<vector>& U_;

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
         * \brief Update thermodynamic fields on boundary face Gauss points
         *
         * \param cellID DG cell index
         */
        virtual void updateBC(const label& cellID) = 0;

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

        /** Return energy model */
        inline const energy& energyModel() const { return energy_(); }

        /** Return gas constant */
        inline scalar R() const { return R_; }

        /** Return stored energy variable name ("e" or "h") */
        inline const word& heName() const { return energy_().he(); }

        /** Return true if he stores sensible internal energy */
        inline bool heIsInternalEnergy() const
        {
            return energy_().heIsInternalEnergy();
        }

        /** Return true if he stores sensible enthalpy */
        inline bool heIsEnthalpy() const
        {
            return energy_().heIsEnthalpy();
        }

        /** Return stored thermodynamic energy field */
        inline const dgField<scalar>& he() const { return he_; }

        /** Return speed-of-sound field */
        inline const dgField<scalar>& a() const { return a_; }

        /** Return enthalpy reconstructed from the stored energy variable */
        tmp<GaussField<scalar>> calcH(const label cellID) const;

        /** Return internal energy reconstructed from the stored energy variable */
        tmp<GaussField<scalar>> calcE(const label cellID) const;

        /** Return pressure field */
        inline dgField<scalar>& p() { return p_; }
        inline const dgField<scalar>& p() const { return p_; }

        /** Return temperature field */
        inline dgField<scalar>& T() { return T_; }
        inline const dgField<scalar>& T() const { return T_; }

        /** Return velocity field */
        inline dgField<vector>& U() { return U_; }
        inline const dgField<vector>& U() const { return U_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
