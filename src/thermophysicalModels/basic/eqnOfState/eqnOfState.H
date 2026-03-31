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
    Foam::eqnOfState

Description
    Base class for Equation of State (single-phase).

SourceFiles
    eqnOfState.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_eqnOfState_H
#define Foam_eqnOfState_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"
#include "dgField.H"
#include "boundaryGaussField.H"
#include "dgGeomMesh.H"
#include "dgMath.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class eqnOfState Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::eqnOfState
 * \brief Base class for single-phase equations of state
 *
 * \details
 * This class defines the interface for single-phase equations of state (EOS)
 * used in DGFoam. An EOS provides thermodynamic closure relations between
 * pressure, temperature, and density, and supplies classification flags
 * that allow thermo models to verify physical compatibility.
 *
 * \par EOS classification and compatibility
 * To ensure consistency between EOS models and thermoLaw formulations,
 * each EOS advertises a set of capability and classification flags.
 * Higher-level thermo models must check these flags in their constructors
 * to ensure that the selected EOS is physically compatible with their
 * assumptions.
 *
 * EOS models are grouped into the following categories:
 *
 * \par Ideal-gas–type EOS
 * - Identified by isIdealGas() returning true.
 * - Pressure relation of the form: \f$ p = \rho R T \f$.
 * - Includes classic ideal gas and thermal-perfect gas variants.
 * - Thermodynamic properties (Cp, Cv, e, h) are provided by thermoLaw.
 * - Compatible with most thermoLaw models.
 *
 * \par Incompressible EOS
 * - Identified by isIncompressible() returning true.
 * - Density is constant or independent of pressure and temperature.
 * - Pressure does not participate in thermodynamic relations.
 * - Only compatible with thermo models designed for liquids.
 *
 * \par Real-gas EOS
 * - Identified by isRealGas() returning true.
 * - Non-ideal equations of state (e.g. Peng–Robinson, Redlich–Kwong,
 *   stiffened gas, Van der Waals).
 * - Require thermo models derived for real-gas properties.
 *
 * \par Kinetic-theory EOS
 * - Identified by isThermalPerfectGas() returning true.
 * - EOS derived from kinetic theory (Boltzmann / ES-BGK).
 * - Pressure may originate from stress rather than \f$ p = \rho R T \f$.
 * - Used for rarefied or high-Knudsen-number flows.
 *
 * By default, all classification flags return false in this base class.
 * Derived EOS models override the appropriate flags to identify their
 * physical category.
 *
 * \ingroup grpDGThermo
 */
class eqnOfState
{
protected:

    /** Model name (used for logging and diagnostics) */
    word name_;

    /** Coefficients dictionary (owned, immutable snapshot) */
    const dictionary dict_;

    /** Reference to DG geometric mesh (non-owning) */
    const dgGeomMesh& mesh_;

public:

    //- Runtime type information
    TypeName("eqnOfState");

    //- Runtime selection table
    declareRunTimeSelectionTable
    (
        autoPtr, eqnOfState, dictionary,
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        ),
        (name, dict, mesh)
    );

    /**
     * \brief Construct from model name and coefficient dictionary
     *
     * \param name Model name
     * \param dict Coefficient dictionary
     * \param mesh DG geometric mesh
     */
    eqnOfState
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** Virtual destructor */
    virtual ~eqnOfState() = default;

    /**
     * \brief Runtime selector
     *
     * Construct an equation-of-state model from a dictionary entry.
     *
     * \param name Model name
     * \param dict Coefficient dictionary
     * \param mesh DG geometric mesh
     * \return Auto pointer to selected EOS model
     */
    static autoPtr<eqnOfState> New
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //====================================================================//
    //                    General EOS interface
    //====================================================================//

    /**
     * \brief Return gas constant R
     *
     * Returns zero for incompressible EOS models.
     */
    virtual scalar R() const = 0;


    //====================================================================//
    //              Fundamental thermodynamic relations (pointwise)
    //
    // Derived EOS models implement scalar pointwise relations only.
    // Field-level overloads below remain as compatibility wrappers that
    // loop over Gauss points and dispatch to these scalar functions.
    //====================================================================//

    /**
     * \brief Compute density from pressure and temperature at one point
     */
    virtual scalar calcRhoFromPT
    (
        const scalar p,
        const scalar T
    ) const
    {
        FatalErrorInFunction
            << "Pointwise calcRhoFromPT() not supported by this EOS model"
            << exit(FatalError);

        return Zero;
    }

    /**
     * \brief Compute pressure from density and temperature at one point
     */
    virtual scalar calcPFromRhoT
    (
        const scalar rho,
        const scalar T
    ) const
    {
        FatalErrorInFunction
            << "Pointwise calcPFromRhoT() not supported by this EOS model"
            << exit(FatalError);

        return Zero;
    }

    /**
     * \brief Compute temperature from pressure and density at one point
     */
    virtual scalar calcTFromPRho
    (
        const scalar p,
        const scalar rho
    ) const
    {
        FatalErrorInFunction
            << "Pointwise calcTFromPRho() not supported by this EOS model"
            << exit(FatalError);

        return Zero;
    }

    //====================================================================//
    //           Fundamental thermodynamic relations (Gauss fields)
    //
    // These wrappers preserve the existing field-based API while the
    // concrete model implementations stay fully pointwise.
    //====================================================================//

    /**
     * \brief Compute density from pressure and temperature
     *
     * Computes \f$ \rho = \rho(p, T) \f$ at Gauss points of a cell.
     *
     * \param cellI Cell index
     * \param p     Pressure field
     * \param T     Temperature field
     * \param rho   Density field to be written
     */
    virtual void calcRhoFromPT
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& T,
        GaussField<scalar>& rho
    ) const;

    virtual void calcRhoFromPT
    (
        const boundaryGaussField<scalar>& p,
        const boundaryGaussField<scalar>& T,
        boundaryGaussField<scalar>& rho
    ) const;

    /**
     * \brief Compute pressure from density and temperature
     *
     * Computes \f$ p = p(\rho, T) \f$ at Gauss points of a cell.
     */
    virtual void calcPFromRhoT
    (
        const label cellI,
        const GaussField<scalar>& rho,
        const GaussField<scalar>& T,
        GaussField<scalar>& p
    ) const;

    virtual void calcPFromRhoT
    (
        const boundaryGaussField<scalar>& rho,
        const boundaryGaussField<scalar>& T,
        boundaryGaussField<scalar>& p
    ) const;

    /**
     * \brief Compute temperature from pressure and density
     *
     * Computes \f$ T = T(p, \rho) \f$ at Gauss points of a cell.
     */
    virtual void calcTFromPRho
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& rho,
        GaussField<scalar>& T
    ) const;

    virtual void calcTFromPRho
    (
        const boundaryGaussField<scalar>& p,
        const boundaryGaussField<scalar>& rho,
        boundaryGaussField<scalar>& T
    ) const;

    //====================================================================//
    //                    EOS capability flags
    //====================================================================//

    /** Return true if rho(p,T) is supported */
    virtual bool canCalcRhoFromPT() const { return false; }

    /** Return true if p(rho,T) is supported */
    virtual bool canCalcPFromRhoT() const { return false; }

    /** Return true if T(p,rho) is supported */
    virtual bool canCalcTFromPRho() const { return false; }

    // EOS classification flags

        /** Ideal-gas–type EOS */
        virtual bool isIdealGas() const { return false; }

        /** Incompressible EOS */
        virtual bool isIncompressible() const { return false; }

        /** Real-gas EOS */
        virtual bool isRealGas() const { return false; }

        /** Kinetic-theory / thermal-perfect EOS */
        virtual bool isThermalPerfectGas() const { return false; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
