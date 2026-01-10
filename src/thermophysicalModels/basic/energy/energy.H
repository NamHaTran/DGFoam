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

/**
 * \class Foam::energy
 * \brief Base class for energy models in DGFoam
 *
 * \details
 * Energy models define which thermodynamic energy variable is evolved
 * by the DG solver and provide consistent conversions between temperature
 * and energy representations.
 *
 * Two categories of energy formulations are supported:
 *
 * 1. Internal-energy based models
 * - The conserved variable is the sensible internal energy e(T).
 * - Temperature is obtained from e(T) using thermoLaw relations.
 * - Enthalpy is computed as: h = e + p/rho.
 * - No formation enthalpy is included.
 * - Suitable for compressible flows, shock-capturing, DG solvers,
 *   and real-gas thermodynamics.
 *
 * 2. Enthalpy-based models
 * - The conserved variable is the sensible enthalpy h(T).
 * - Temperature is obtained from h(T).
 * - Internal energy is computed as: e = h - p/rho.
 * - Required when transport or chemistry formulations are expressed
 *   naturally in terms of enthalpy.
 *
 * This base class provides a common interface for all energy models
 * and supports runtime selection via dictionaries.
 *
 * \ingroup grpDGThermo
 */
class energy
{
protected:

    /** Model name (used for logging and diagnostics) */
    word name_;

    /** Coefficient dictionary snapshot */
    const dictionary coeffs_;

    /** Reference to DG geometric mesh (non-owning) */
    const dgGeomMesh& mesh_;

    /** Reference to thermoLaw (non-owning) */
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

    /**
     * \brief Construct from name, dictionary, mesh, and thermoLaw
     *
     * \param name   Model name
     * \param dict   Coefficient dictionary
     * \param mesh   DG geometric mesh
     * \param thermo Thermodynamic law
     */
    energy
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    /** Virtual destructor */
    virtual ~energy() = default;

    /**
     * \brief Runtime selector
     *
     * \param name   Model name
     * \param dict   Coefficient dictionary
     * \param mesh   DG geometric mesh
     * \param thermo Thermodynamic law
     * \return Auto pointer to selected energy model
     */
    static autoPtr<energy> New
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    // * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

    /**
     * \brief Compute enthalpy from temperature
     *
     * The computed enthalpy is written directly into the output field
     * for the specified cell.
     *
     * \param cellI Cell index
     * \param T     Temperature field at Gauss points
     * \param h     Enthalpy field to be written
     */
    virtual void calcEnthalpy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& h
    ) const = 0;

    /**
     * \brief Compute energy variable (e or h depending on model)
     *
     * \param cellI Cell index
     * \param T     Temperature field at Gauss points
     * \param E     Energy field to be written
     */
    virtual void calcEnergy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& E
    ) const = 0;

    /**
     * \brief Compute temperature from energy (energy inversion)
     *
     * \param cellI Cell index
     * \param E     Energy field at Gauss points
     * \param T     Temperature field to be written
     */
    virtual void calcTfromEnergy
    (
        const label cellI,
        const GaussField<scalar>& E,
        GaussField<scalar>& T
    ) const = 0;

    // Energy model identification

        /** Return true if the model is internal-energy based */
        virtual bool energyInternal() const { return false; }

        /** Return true if the model is enthalpy-based */
        virtual bool energyEnthalpy() const { return false; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
