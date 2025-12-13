/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DGFoam: Discontinuous Galerkin CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | GPU-friendly CFD solver framework
     \\/     M anipulation  |
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgFluxSolverManager.H"
#include "IOstreams.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dgFluxSolverManager::dgFluxSolverManager
(
    const dictionary& dgSchemesDict,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    fluxSolverList_(),
    termList_(),
    defaultFluxSolver_(nullptr),
    // Cache the whole block for optional later use
    fluxSchemesDict_(dgSchemesDict.subDict("fluxSchemes"))
{
    // 0) Locate sub-dicts: fluxSolvers & fluxSolversCoeffs
    if (!fluxSchemesDict_.found("fluxSolvers"))
    {
        FatalIOErrorInFunction(dgSchemesDict)
            << "Missing sub-dictionary 'fluxSchemes/fluxSolvers'." << nl
            << exit(FatalIOError);
    }

    const dictionary& solversDict = fluxSchemesDict_.subDict("fluxSolvers");

    const dictionary* coeffsRootPtr = nullptr;
    if (fluxSchemesDict_.found("fluxSolversCoeffs"))
    {
        coeffsRootPtr = &fluxSchemesDict_.subDict("fluxSolversCoeffs");
    }
    else
    {
        FatalIOErrorInFunction(dgSchemesDict)
            << "Missing sub-dictionary 'fluxSchemes/fluxSolversCoeffs'." << nl
            << exit(FatalIOError);
    }

    // 1) Default solver (optional)
    if (solversDict.found("default"))
    {
        const word defaultScheme(solversDict.get<word>("default"));

        dictionary defSolverDict; // empty by default
        if (coeffsRootPtr)
        {
            const word coeffKey(defaultScheme + "Coeffs");
            if (coeffsRootPtr->found(coeffKey))
            {
                defSolverDict = coeffsRootPtr->subDict(coeffKey); // copy
            }
            else
            {
                FatalIOErrorInFunction(dgSchemesDict)
                    << "Missing sub-dictionary 'fluxSchemes/fluxSolversCoeffs/" << coeffKey << "'." << nl
                    << exit(FatalIOError);
            }
        }

        defaultFluxSolver_.reset
        (
            dgFluxSolver::New("default", defaultScheme, defSolverDict, mesh).ptr()
        );
    }
    // else: keep defaultFluxSolver_ = nullptr

    // 2) Collect term entries: skip "default"
    DynamicList<word> terms;
    forAllConstIter(dictionary, solversDict, iter)
    {
        const word key = iter().keyword();

        if (key == "FoamFile") continue; // ignore header if present
        if (key == "default")  continue; // skip default

        // Only accept simple "term scheme;" entries (not nested dicts)
        if (solversDict.isDict(key)) continue;

        terms.append(key);
    }
    termList_.transfer(terms);

    // 3) Instantiate solvers for each configured term
    fluxSolverList_.setSize(termList_.size());

    forAll(termList_, i)
    {
        const word& term   = termList_[i];
        const word  scheme = solversDict.get<word>(term);

        // dict cho solver = block <SchemeName>Coeffs (if any)
        dictionary sDict; // empty by default
        if (coeffsRootPtr)
        {
            const word coeffKey(scheme + "Coeffs");
            if (coeffsRootPtr->found(coeffKey))
            {
                sDict = coeffsRootPtr->subDict(coeffKey); // copy
            }
            else
            {
                FatalIOErrorInFunction(dgSchemesDict)
                    << "Missing sub-dictionary 'fluxSchemes/fluxSolversCoeffs/" << coeffKey << "'." << nl
                    << exit(FatalIOError);
            }
        }

        autoPtr<dgFluxSolver> slv =
            dgFluxSolver::New(term, scheme, sDict, mesh);

        fluxSolverList_.set(i, slv.ptr());
    }
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

inline bool dgFluxSolverManager::has(const word& term) const
{
    forAll(termList_, i)
    {
        if (termList_[i] == term) return true;
    }
    return false;
}


const dgFluxSolver& dgFluxSolverManager::solver(const word& term) const
{
    forAll(termList_, i)
    {
        if (termList_[i] == term)
        {
            return fluxSolverList_[i];
        }
    }

    if (defaultFluxSolver_.valid())
    {
        return defaultFluxSolver_();
    }

    FatalErrorInFunction
        << "No flux scheme configured for term \"" << term << "\" and no "
        << "default scheme provided.\n"
        << "Available terms: " << termList_ << nl
        << exit(FatalError);

    return fluxSolverList_[0]; // Unreachable
}


dgFluxSolver& dgFluxSolverManager::solver(const word& term)
{
    forAll(termList_, i)
    {
        if (termList_[i] == term)
        {
            return fluxSolverList_[i];
        }
    }

    if (defaultFluxSolver_.valid())
    {
        return defaultFluxSolver_();
    }

    FatalErrorInFunction
        << "No flux scheme configured for term \"" << term << "\" and no "
        << "default scheme provided.\n"
        << "Provided flux schemes: " << termList_ << nl
        << exit(FatalError);

    return fluxSolverList_[0]; // Unreachable
}


void dgFluxSolverManager::listTerms(Ostream& os) const
{
    os  << "dgFluxSolverManager: configured terms = "
        << termList_.size() << nl;

    forAll(termList_, i)
    {
        const word& term = termList_[i];
        // Scheme name is stored in each solver's dict under key "solver"
        word scheme("unknown");
        const dictionary& sd = fluxSolverList_[i].dict();
        if (sd.found("solver")) scheme = sd.get<word>("solver");

        os << "  - " << term << " : " << scheme << nl;
    }

    if (defaultFluxSolver_.valid())
    {
        word scheme("unknown");
        const dictionary& sd = defaultFluxSolver_().dict();
        if (sd.found("solver")) scheme = sd.get<word>("solver");

        os << "  default : " << scheme << nl;
    }
    else
    {
        os << "  default : (none)" << nl;
    }
}

} // End namespace Foam

// ************************************************************************* //

