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
    defaultConvectiveFluxSolver_(nullptr),
    defaultDiffusiveFluxSolver_(nullptr),
    // Cache the whole block so each solver can read its own coeff sub-dict
    fluxSchemesDict_(dgSchemesDict.subDict("fluxSchemes"))
{
    // 0) Locate sub-dict: fluxSolvers
    if (!fluxSchemesDict_.found("fluxSolvers"))
    {
        FatalIOErrorInFunction(dgSchemesDict)
            << "Missing sub-dictionary 'fluxSchemes/fluxSolvers'." << nl
            << exit(FatalIOError);
    }

    const dictionary& solversDict = fluxSchemesDict_.subDict("fluxSolvers");

    if (solversDict.found("default"))
    {
        FatalIOErrorInFunction(solversDict)
            << "The fluxSchemes/fluxSolvers keyword 'default' is no longer "
            << "supported.\n"
            << "Use 'defaultConvective' and/or 'defaultDiffusive' instead."
            << nl << exit(FatalIOError);
    }

    // 1) Default solvers (optional)
    if (solversDict.found("defaultConvective"))
    {
        const word defaultScheme(solversDict.get<word>("defaultConvective"));

        defaultConvectiveFluxSolver_.reset
        (
            dgFluxSolver::New
            (
                "defaultConvective",
                defaultScheme,
                fluxSchemesDict_,
                mesh
            ).ptr()
        );
    }

    if (solversDict.found("defaultDiffusive"))
    {
        const word defaultScheme(solversDict.get<word>("defaultDiffusive"));

        defaultDiffusiveFluxSolver_.reset
        (
            dgFluxSolver::New
            (
                "defaultDiffusive",
                defaultScheme,
                fluxSchemesDict_,
                mesh
            ).ptr()
        );
    }

    // 2) Collect term entries: skip defaults
    DynamicList<word> terms;
    forAllConstIter(dictionary, solversDict, iter)
    {
        const word key = iter().keyword();

        if (key == "FoamFile")          continue; // ignore header if present
        if (key == "defaultConvective") continue; // skip default
        if (key == "defaultDiffusive")  continue; // skip default

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

        autoPtr<dgFluxSolver> slv =
            dgFluxSolver::New(term, scheme, fluxSchemesDict_, mesh);

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


const dgFluxSolver& dgFluxSolverManager::solver
(
    const word& term,
    dgFluxSolver::fluxType fType
) const
{
    const dgFluxSolver* slvPtr = nullptr;

    forAll(termList_, i)
    {
        if (termList_[i] == term)
        {
            slvPtr = &fluxSolverList_[i];
            break;
        }
    }

    if (!slvPtr)
    {
        if
        (
            fType == dgFluxSolver::fluxType::convective
         && defaultConvectiveFluxSolver_.valid()
        )
        {
            slvPtr = &defaultConvectiveFluxSolver_();
        }
        else if
        (
            fType == dgFluxSolver::fluxType::diffusive
         && defaultDiffusiveFluxSolver_.valid()
        )
        {
            slvPtr = &defaultDiffusiveFluxSolver_();
        }
    }

    if (!slvPtr)
    {
        FatalErrorInFunction
            << "No flux scheme configured for term \"" << term << "\" and no "
            << fType << " default scheme provided.\n"
            << "Available terms: " << termList_ << nl
            << exit(FatalError);
    }

    const dgFluxSolver& slv = *slvPtr;

    if
    (
        slv.getFluxType() != fType
     && slv.getFluxType() != dgFluxSolver::fluxType::both
    )
    {
        FatalErrorInFunction
            << "Flux solver \"" << slv.name() << "\" selected for term \""
            << term << "\" has flux type " << slv.getFluxType()
            << ", but the PDE term requested " << fType << ".\n"
            << "Check the fluxSchemes/fluxSolvers entry for this term."
            << nl << exit(FatalError);
    }

    return slv;
}


dgFluxSolver& dgFluxSolverManager::solver
(
    const word& term,
    dgFluxSolver::fluxType fType
)
{
    dgFluxSolver* slvPtr = nullptr;

    forAll(termList_, i)
    {
        if (termList_[i] == term)
        {
            slvPtr = &fluxSolverList_[i];
            break;
        }
    }

    if (!slvPtr)
    {
        if
        (
            fType == dgFluxSolver::fluxType::convective
         && defaultConvectiveFluxSolver_.valid()
        )
        {
            slvPtr = &defaultConvectiveFluxSolver_();
        }
        else if
        (
            fType == dgFluxSolver::fluxType::diffusive
         && defaultDiffusiveFluxSolver_.valid()
        )
        {
            slvPtr = &defaultDiffusiveFluxSolver_();
        }
    }

    if (!slvPtr)
    {
        FatalErrorInFunction
            << "No flux scheme configured for term \"" << term << "\" and no "
            << fType << " default scheme provided.\n"
            << "Provided flux schemes: " << termList_ << nl
            << exit(FatalError);
    }

    dgFluxSolver& slv = *slvPtr;

    if
    (
        slv.getFluxType() != fType
     && slv.getFluxType() != dgFluxSolver::fluxType::both
    )
    {
        FatalErrorInFunction
            << "Flux solver \"" << slv.name() << "\" selected for term \""
            << term << "\" has flux type " << slv.getFluxType()
            << ", but the PDE term requested " << fType << ".\n"
            << "Check the fluxSchemes/fluxSolvers entry for this term."
            << nl << exit(FatalError);
    }

    return slv;
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

    if (defaultConvectiveFluxSolver_.valid())
    {
        word scheme("unknown");
        const dictionary& sd = defaultConvectiveFluxSolver_().dict();
        if (sd.found("solver")) scheme = sd.get<word>("solver");

        os << "  defaultConvective : " << scheme << nl;
    }
    else
    {
        os << "  defaultConvective : (none)" << nl;
    }

    if (defaultDiffusiveFluxSolver_.valid())
    {
        word scheme("unknown");
        const dictionary& sd = defaultDiffusiveFluxSolver_().dict();
        if (sd.found("solver")) scheme = sd.get<word>("solver");

        os << "  defaultDiffusive : " << scheme << nl;
    }
    else
    {
        os << "  defaultDiffusive : (none)" << nl;
    }
}

} // End namespace Foam

// ************************************************************************* //
