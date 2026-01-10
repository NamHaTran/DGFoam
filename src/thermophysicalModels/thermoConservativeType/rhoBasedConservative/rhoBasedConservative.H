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
 * The reconstructed quantities include:
 * \li Temperature \f$T\f$
 * \li Pressure \f$p\f$
 * \li Internal energy \f$e\f$
 * \li Enthalpy \f$h\f$
 * \li Heat capacities \f$C_p\f$, \f$C_v\f$
 * \li Heat-capacity ratio \f$\gamma\f$
 * \li Speed of sound \f$a\f$
 * \li Transport properties \f$\mu\f$, \f$\kappa\f$, \f$\mathrm{Pr}\f$
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
};

} // End namespace Foam

#endif

// ************************************************************************* //
