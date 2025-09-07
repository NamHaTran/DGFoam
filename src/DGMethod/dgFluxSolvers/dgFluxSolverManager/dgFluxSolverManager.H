/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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
    Foam::dgFluxSolverManager

Description
    Manager for flux solvers, mapping equation terms (e.g. massConvectiveFlux)
    to concrete dgFluxSolver instances. Provides a default solver when a term
    has no explicit scheme in the dictionary.

    Dictionary format (sub-dict "fluxSchemes"):
        default                 pureCentral;
        massConvectiveFlux      HLLE;
        momentumConvectiveFlux  LF;
        HLLECoeff
        {
            speedEstimate roeEinfeldt;
        }

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
#include "dgFluxSolver.H"
#include "dgThermo.H"

namespace Foam
{

// Forward declare Ostream to avoid heavy includes here
class Ostream;

/*---------------------------------------------------------------------------*\
                      Class dgFluxSolverManager Declaration
\*---------------------------------------------------------------------------*/

class dgFluxSolverManager
{
private:

    // Attributes

        //- List of flux solver instances (one per configured term)
        PtrList<dgFluxSolver> fluxSolverList_;

        //- Parallel list of term names corresponding to fluxSolverList_ items
        List<word> termList_;

        //- Default flux solver used when a term is not configured (may be null)
        autoPtr<dgFluxSolver> defaultFluxSolver_;

        //- Cached schemes dictionary (fluxSchemes)
        dictionary fluxSchemesDict_;

public:

    // Constructors

        //- Construct from dgSchemes dictionary and thermo (build solvers)
        dgFluxSolverManager
        (
            const dictionary& dgSchemesDict,
            const dgThermo& thermo
        );

        //- Disallow default construction
        dgFluxSolverManager() = delete;

        //- Destructor
        ~dgFluxSolverManager() = default;

    // Member Functions

        //- Set context for all managed solvers (including default if any)
        void setContext(const fieldsContext& ctx);

        //- Number of configured (term-specific) solvers
        inline label size() const { return fluxSolverList_.size(); }

        //- Check if a term has an explicitly configured solver (no fallback)
        inline bool has(const word& term) const;

        //- Access solver by term name; falls back to default or throws on missing
        const dgFluxSolver& solver(const word& term) const;

        //- Non-const access variant
        dgFluxSolver& solver(const word& term);

        //- Print configured terms and their scheme names (and default) to stream
        void listTerms(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //