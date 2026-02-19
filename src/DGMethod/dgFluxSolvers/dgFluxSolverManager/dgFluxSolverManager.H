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
    Foam::dgFluxSolverManager

Description
    Manager for DG flux solvers selected per equation term.

SourceFiles
    dgFluxSolverManager.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFluxSolverManager_H
#define Foam_dgFluxSolverManager_H

#include "PtrList.H"
#include "List.H"
#include "word.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "dgGeomMesh.H"
#include "dgFluxSolver.H"

namespace Foam
{

// Forward declare Ostream to avoid heavy includes here
class Ostream;

/*---------------------------------------------------------------------------*\
                      Class dgFluxSolverManager Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgFluxSolverManager
 * \brief Manager for flux solvers selected per equation term
 *
 * \details
 * This class constructs and manages dgFluxSolver instances based on
 * the \c fluxSchemes dictionary.
 *
 * Each equation term (e.g. massConvTermFlux, momentumConvTermFlux)
 * can be mapped to a specific numerical flux solver. If a term is
 * not explicitly defined, a default solver is used.
 *
 * The manager supports:
 * - Term-specific flux solver selection
 * - A global default solver fallback
 * - Solver-specific coefficient sub-dictionaries
 *
 * \par Example dictionary
 * \verbatim
 * fluxSchemes
 * {
 *     fluxSolvers
 *     {
 *         default                 central;
 *
 *         massConvTermFlux        HLLC;
 *         momentumConvTermFlux    LaxFriedrichs;
 *
 *         // Other flux terms such as:
 *         // - auxiliaryFlux
 *         // - massDiffusiveFlux
 *         // - momentumDiffusiveFlux
 *         // - energyConvectiveFlux
 *         // - energyDiffusiveFlux
 *         // fall back to the default solver if not specified.
 *     }
 *
 *     fluxSolversCoeffs
 *     {
 *         HLLCCoeffs
 *         {
 *             speedEstimate roeEinfeldt;
 *         }
 *
 *         LaxFriedrichsCoeffs
 *         {
 *             scaleByMach true;
 *         }
 *
 *         centralCoeffs
 *         {
 *             // no parameters
 *         }
 *     }
 * }
 * \endverbatim
 *
 * Solver instances are created lazily using runtime selection tables
 * and are reused across time steps. The reset() method propagates to
 * all managed solvers.
 */
class dgFluxSolverManager
{
private:

    // Attributes

        /** Reference to DG geometric mesh */
        const dgGeomMesh& mesh_;

        /** Flux solver instances (one per configured term) */
        PtrList<dgFluxSolver> fluxSolverList_;

        /** Term names corresponding to fluxSolverList_ ordering */
        List<word> termList_;

        /** Default flux solver used when a term is not configured (may be null) */
        autoPtr<dgFluxSolver> defaultFluxSolver_;

        /** Cached \c fluxSchemes sub-dictionary */
        dictionary fluxSchemesDict_;

public:

    // Constructors

        /**
         * \brief Construct from dgSchemes dictionary and DG mesh
         *
         * Parses \c fluxSchemes and constructs term-specific solvers and
         * an optional default solver.
         *
         * \param dgSchemesDict dgSchemes dictionary
         * \param mesh          DG geometric mesh
         */
        dgFluxSolverManager
        (
            const dictionary& dgSchemesDict,
            const dgGeomMesh& mesh
        );

        //- Disallow default construction
        dgFluxSolverManager() = delete;

        //- Destructor
        ~dgFluxSolverManager() = default;

    // Member Functions

        /** Return number of term-specific solvers */
        inline label size() const { return fluxSolverList_.size(); }

        /**
         * \brief Check if a term has an explicitly configured solver
         *
         * \param term Equation-term name
         * \return true if configured in \c fluxSchemes
         */
        inline bool has(const word& term) const;

        /**
         * \brief Return solver by term name (const)
         *
         * Falls back to the default solver if available; otherwise triggers
         * an error for missing configuration.
         *
         * \param term Equation-term name
         * \return Const reference to selected solver
         */
        const dgFluxSolver& solver(const word& term) const;

        /**
         * \brief Reset all solvers
         *
         * Calls reset() on all term-specific solvers and the default solver
         * (if present).
         */
        inline void reset()
        {
            forAll(fluxSolverList_, i)
            {
                fluxSolverList_[i].reset();
            }

            if (defaultFluxSolver_.valid())
            {
                defaultFluxSolver_.reset();
            }
        }

        /**
         * \brief Return solver by term name (non-const)
         *
         * \param term Equation-term name
         * \return Reference to selected solver
         */
        dgFluxSolver& solver(const word& term);

        /**
         * \brief Print configured terms and their solver types
         *
         * Writes the mapping of term names to solver types, including the
         * default solver if configured.
         *
         * \param os Output stream
         */
        void listTerms(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
