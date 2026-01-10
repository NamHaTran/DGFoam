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
    Foam::thermoLaw

Description
    Base class for thermodynamic laws used in DG thermophysical models.

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

/**
 * \class Foam::thermoLaw
 * \brief Base class for thermodynamic laws in DG solvers
 *
 * \details
 * thermoLaw defines the thermodynamic closure relations required by DG
 * thermophysical models. It provides a unified interface for:
 *
 * \li Converting between temperature and energy forms \f$(e, h)\f$
 * \li Computing heat capacities \f$C_p(T)\f$, \f$C_v(T)\f$
 * \li Computing heat-capacity ratio \f$\gamma(T)\f$
 * \li Computing speed of sound \f$a(T)\f$
 * \li Evaluating internal energy and enthalpy at Gauss points
 *
 * All operations are performed locally at Gauss points of each DG cell,
 * consistent with high-order discontinuous Galerkin formulations.
 *
 * \par Thermodynamic model categories
 *
 * \par Perfect-gas thermodynamics
 * \li \f$C_p(T)\f$ and \f$C_v(T)\f$ depend only on temperature
 * \li No formation enthalpy is included
 * \li \f$h(T) = \int C_p(T)\,\mathrm{d}T\f$
 * \li \f$e(T) = \int C_v(T)\,\mathrm{d}T\f$
 * \li Compatible with ideal-gas EOS: \f$p = \rho R T\f$
 * \li Used in compressible and high-temperature gas flows
 *
 * \par Absolute-enthalpy thermodynamics
 * \li Enthalpy includes chemical formation energy
 * \li \f$h(T) = h_f^0 + \int C_p(T)\,\mathrm{d}T\f$
 * \li Required for combustion and reacting flows
 * \li Works with ideal-gas or selected real-gas EOS
 *
 * \par Real-gas thermodynamics
 * \li \f$C_p\f$, \f$C_v\f$, \f$e\f$, and \f$h\f$ depend on both temperature
 *     and density
 * \li Must be used together with real-gas EOS
 * \li Required for high-pressure and supercritical fluids
 *
 * \par Kinetic / rarefied-gas thermodynamics
 * \li Based on kinetic theory and molecular degrees of freedom
 * \li \f$e = \frac{f}{2} R T\f$
 * \li \f$h = \frac{f+2}{2} R T\f$
 * \li \f$C_p\f$ and \f$C_v\f$ follow kinetic DOF relations
 * \li Compatible with thermal-perfect-gas EOS
 * \li Used in rarefied gas and hypersonic applications
 */
class thermoLaw
{
protected:

    /** Model name (used for logging and diagnostics) */
    word name_;

    /** Coefficient dictionary (immutable after construction) */
    const dictionary dict_;

    /** Reference to DG geometric mesh (non-owning) */
    const dgGeomMesh& mesh_;

    /** Reference to equation-of-state model */
    const eqnOfState& eos_;

public:

    //- Runtime type information
    TypeName("thermoLaw");

    //- Runtime selection table
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

    /**
     * \brief Construct from model name and coefficient dictionary
     *
     * \param name Model name
     * \param dict Coefficient dictionary
     * \param mesh DG geometric mesh
     * \param eos  Equation-of-state model
     */
    thermoLaw
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const eqnOfState& eos
    );

    /** Virtual destructor */
    virtual ~thermoLaw() = default;

    /**
     * \brief Runtime selector
     *
     * \param name Model name
     * \param dict Coefficient dictionary
     * \param mesh DG geometric mesh
     * \param eos  Equation-of-state model
     * \return Auto pointer to selected thermoLaw
     */
    static autoPtr<thermoLaw> New
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const eqnOfState& eos
    );

    //====================================================================//
    //                   General thermodynamic interface
    //====================================================================//

    /**
     * \brief Compute temperature from energy
     *
     * Computes temperature $T$ from internal energy $e(T)$ or enthalpy $h(T)$,
     * depending on the specific thermoLaw model.
     *
     * \param cellI  Cell index
     * \param energy Energy field ($e$ or $h$)
     * \param T      Temperature field to be written
     */
    virtual void calcT
    (
        const label cellI,
        const GaussField<scalar>& energy,
        GaussField<scalar>& T
    ) const = 0;

    /**
     * \brief Compute heat capacity at constant pressure
     *
     * \param cellI Cell index
     * \param T     Temperature field
     * \param Cp    $C_p(T)$ field
     */
    virtual void calcCp
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& Cp
    ) const = 0;

    /**
     * \brief Compute heat capacity at constant volume
     *
     * \param cellI Cell index
     * \param T     Temperature field
     * \param Cv    $C_v(T)$ field
     */
    virtual void calcCv
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& Cv
    ) const = 0;

    /**
     * \brief Compute heat-capacity ratio
     *
     * Computes $\gamma = C_p / C_v$.
     */
    virtual void calcGamma
    (
        const label cellI,
        const GaussField<scalar>& Cp,
        const GaussField<scalar>& Cv,
        GaussField<scalar>& gamma
    ) const = 0;

    /**
     * \brief Compute internal energy
     *
     * Computes internal energy $e(T)$.
     */
    virtual void calcInternalE
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& e
    ) const = 0;

    /**
     * \brief Compute enthalpy
     *
     * Computes enthalpy $h(T)$.
     */
    virtual void calcH
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& h
    ) const = 0;

    /**
     * \brief Compute speed of sound
     *
     * Computes speed of sound $a = \sqrt{\gamma R T}$ or the
     * corresponding real-gas / kinetic-theory expression.
     */
    virtual void calcSpeedOfSound
    (
        const label cellI,
        const GaussField<scalar>& T,
        const GaussField<scalar>& gamma,
        GaussField<scalar>& a
    ) const = 0;

    //====================================================================//
    //                   ThermoLaw classification flags
    //====================================================================//

        /** Perfect-gas thermo (temperature-dependent Cp, no formation energy) */
        virtual bool isPerfectGasThermo() const { return false; }

        /** Absolute-enthalpy thermo (includes formation enthalpy) */
        virtual bool isAbsoluteEnthalpyThermo() const { return false; }

        /** Real-gas thermo */
        virtual bool isRealGasThermo() const { return false; }

        /** Kinetic / rarefied-gas thermo */
        virtual bool isKineticThermo() const { return false; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
