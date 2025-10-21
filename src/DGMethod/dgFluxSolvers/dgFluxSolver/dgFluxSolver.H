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
    Foam::dgFluxSolver

Description
    Base class for numerical flux solvers that operate on vector-valued
    fluxes. The interface accepts Plus/Minus flux vectors (Cartesian)
    and returns the numerical flux vector.

    Notes
    - New convention:
        * Plus  = neighbour cell (N-side)
        * Minus = owner cell (P-side)
    - This class only works with Foam::vector (no tensor support).
    - A fieldsContext pointer may be set for field lookup/registration.

    Flux solvers's category:
    
    dgFluxSolvers/
    ├── exactRiemannFluxSolver/
    │    ├── Godunov.C
    ├── approximateRiemannFluxSolver/
         ├── centralSchemes/
         │    ├── LaxFriedrichs.C
         │    └── Rusanov.C
         ├── contactResolvingSchemes/
         │    ├── HLLC.C
         │    ├── Roe.C
         └── fluxSplittingSchemes/
              ├── AUSM.C
              ├── VanLeer.C
              └── StegerWarming.C


SourceFiles
    dgFluxSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFluxSolver_H
#define Foam_dgFluxSolver_H

#include "word.H"
#include "dictionary.H"
#include "vector.H"
#include "autoPtr.H"
#include "dgThermo.H"   // Thermo interface for fluid properties

namespace Foam
{

// Forward declaration
class fieldsContext;

/*---------------------------------------------------------------------------*\
                           Class dgFluxSolver Declaration
\*---------------------------------------------------------------------------*/

class dgFluxSolver
{
public:
    // Enumeration for identifying the PDE being solved
    enum class equationType
    {
        massTransport,
        momentumTransport,
        energyTransport,
        scalarTransport
    };

protected:

    // Attributes

        //- Name of flux solver (instance name for logging/debug)
        word name_;

        //- Configuration sub-dictionary for this solver
        dictionary dict_;

        //- Pointer to fields context (for lookup/registration); may be null
        const fieldsContext* ctxPtr_ = nullptr;

        //- Reference to thermo object for fluid properties (non-owning)
        const dgThermo& thermo_;

        // Internal equation type used for interpreting UL/UR
        equationType eqnType_;

public:

    //- Runtime type info
    TypeName("dgFluxSolver");

    //- Runtime selection support
    declareRunTimeSelectionTable
    (
        autoPtr,
        dgFluxSolver,
        dictionary,
        (const word& name, const dictionary& dict, const dgThermo& thermo),
        (name, dict, thermo)
    );

    // Constructors

        //- Construct from name and dictionary
        dgFluxSolver(const word& name, const dictionary& dict, const dgThermo& thermo);

        //- Disallow default construction
        dgFluxSolver() = delete;

        //- Virtual destructor
        virtual ~dgFluxSolver() = default;

    // Factory

        //- Select and construct from dictionary entry "solver"
        static autoPtr<dgFluxSolver> New
        (
            const word& name,
            const word& fluxSolverType,
            const dictionary& dict,
            const dgThermo& thermo
        );

    // Context

        //- Set the fields context (no ownership)
        inline void setContext(const fieldsContext& ctx) { ctxPtr_ = &ctx; }

        //- Get the fields context pointer (may be null)
        inline const fieldsContext* context() const { return ctxPtr_; }

    // Access

        //- Return solver name
        inline const word& name() const { return name_; }

        //- Return solver dictionary
        inline const dictionary& dict() const { return dict_; }

    // Core API

        //- Compute numerical flux from Plus/Minus physical flux vectors
        //  Input:
        //    - gaussID : local Gauss point ID on face
        //    - FPlus   : flux vector from neighbour cell (N-side)
        //    - FMinus  : flux vector from owner cell (P-side)
        //    - UR      : conservative variable scalar from neighbour cell
        //    - UL      : conservative variable scalar from owner cell
        //    - n       : face unit normal (points from P to N)
        //  Output:
        //    - flux    : numerical flux vector (Cartesian)
        virtual void computeFlux
        (
            const label gaussID,
            const vector& FPlus,
            const vector& FMinus,
            const scalar UR,
		    const scalar UL,
            const vector& n,
            vector& flux
        ) const = 0;

    // Utility
        // - Construct orthonormal basis (n, t1, t2) from normal n
        //   (t1, t2 are tangential directions)
        void makeONB(const vector& n, vector& t1, vector& t2) const;

        // - Reading method
        virtual void read(const dictionary& dict) = 0;

        // Equation type Setter and getter
        inline void setEquationType(equationType t) { eqnType_ = t; }
        inline equationType getEquationType() const { return eqnType_; }
};

    // Overload operator<< to print equationType
    inline Ostream& operator<< (Ostream& os, const dgFluxSolver::equationType& type)
    {
        switch (type)
        {
            case dgFluxSolver::equationType::massTransport:
                os << "massTransport"; break;
            case dgFluxSolver::equationType::momentumTransport:
                os << "momentumTransport"; break;
            case dgFluxSolver::equationType::energyTransport:
                os << "energyTransport"; break;
            case dgFluxSolver::equationType::scalarTransport:
                os << "scalarTransport"; break;
            default:
                os << "UnknownEquationType"; break;
        }
        return os;
    }

} // End namespace Foam

#endif

// ************************************************************************* //