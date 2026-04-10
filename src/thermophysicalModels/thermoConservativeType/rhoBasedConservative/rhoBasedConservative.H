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
    Foam::rhoBasedConservative

Description
    Conservative thermodynamics model for DG compressible solvers.

SourceFiles
    rhoBasedConservative.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_rhoBasedConservative_H
#define Foam_rhoBasedConservative_H

#include "dgThermoConservative.H"
#include "eqnOfState.H"
#include "thermoLaw.H"
#include "transportLaw.H"
#include "energy.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class rhoBasedConservative Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::rhoBasedConservative
 * \brief Density-based conservative thermodynamics for DG solvers
 *
 * \details
 * rhoBasedConservative reconstructs all primitive thermodynamic quantities
 * at Gauss points from the conservative variables
 * \f$(\rho, \rho\mathbf{U}, E)\f$.
 *
 * The stored reconstructed quantities include:
 * \li Temperature \f$T\f$
 * \li Pressure \f$p\f$
 * \li Stored energy variable \f$he\f$
 * \li Speed of sound \f$a\f$
 *
 * Additional thermodynamic and transport properties such as
 * \f$C_p\f$, \f$C_v\f$, \f$\gamma\f$, \f$\kappa\f$, \f$\mu\f$, and
 * \f$\mathrm{Pr}\f$ are evaluated on demand through the thermoLaw and
 * transportLaw APIs instead of being stored as dgField data. The meaning of
 * \f$he\f$ is determined by the selected energy model.
 *
 * All operations are performed locally at Gauss points of each DG cell,
 * which is consistent with high-order discontinuous Galerkin discretizations.
 *
 * \par Supported model combinations
 * \li Ideal-gas and thermal-perfect-gas equations of state
 * \li Compatible thermoLaw, transportLaw, and energy models
 *
 * \par Limitations
 * \li Intended for conservative Navier–Stokes formulations
 * \li Not applicable to primitive-variable formulations
 * \li Not applicable to incompressible models
 *
 * \ingroup grpDGThermo
 */
class rhoBasedConservative
:
    public dgThermoConservative
{
private:

        /**
         * \brief Reconstruct a mean-state temperature from conservative mean values.
         *
         * \details
         * This helper converts one cell-mean conservative state
         * \f$(\rho_0,\rho U_0,E_0)\f$ into a fallback mean temperature used
         * when pointwise temperature reconstruction becomes slightly negative.
         */
        scalar calcTMean
        (
            const scalar rho0,
            const vector& rhoU0,
            const scalar E0
        ) const;

public:

    //- Runtime type information
    TypeName("rhoBasedConservative");

    //- Runtime selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoBasedConservative,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        ),
        (name, dict, mesh)
    );

    // Constructors

        /** Disallow default construction */
        rhoBasedConservative() = delete;

        /**
         * \brief Construct from name, dictionary, and DG mesh
         *
         * \param name Model name
         * \param dict Coefficient dictionary
         * \param mesh DG geometric mesh
         */
        rhoBasedConservative
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

    /** Virtual destructor */
    virtual ~rhoBasedConservative() = default;


    // Model setup interface

        /**
         * \brief Construct EOS, thermo, transport, and energy models
         *
         * Creates all sub-models required to reconstruct the thermodynamic
         * state from conservative variables.
         */
        virtual void initModels() override;

        /**
         * \brief Validate model combination
         *
         * Checks physical and mathematical consistency between the selected
         * EOS, thermoLaw, transportLaw, and energy models.
         */
        virtual void validateModels() override;


    // Core API

        /**
         * \brief Update thermodynamic state for a given cell
         *
         * Reconstructs all primitive thermodynamic quantities at Gauss points
         * for the specified cell.
         *
         * \param cellID Cell index
         */
        virtual void update(const label& cellID) override;

        /**
         * \brief Update thermodynamic state on boundary Gauss points
         *
         * \param cellID Cell index
         */
        virtual void updateBC(const label& cellID) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
