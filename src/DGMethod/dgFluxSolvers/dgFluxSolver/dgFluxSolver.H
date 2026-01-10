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
#include "dgGeomMesh.H"
#include "faceGaussField.H"

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

        //- Mesh reference
        const dgGeomMesh& mesh_;

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
        (const word& name, const dictionary& dict, const dgGeomMesh& mesh),
        (name, dict, mesh)
    );

    // Constructors

        //- Construct from name and dictionary
        dgFluxSolver(const word& name, const dictionary& dict, const dgGeomMesh& mesh);

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
            const dgGeomMesh& mesh
        );

    // Access

        //- Return solver name
        inline const word& name() const { return name_; }

        //- Return solver dictionary
        inline const dictionary& dict() const { return dict_; }

    // Core API

        //- Compute numerical flux from Plus/Minus physical flux vectors
        //  Input:
        //    - cellID  : cellID of cell where flux is calculated
        //    - F       : faceGaussField of field F
        //    - U       : faceGaussField of U
        virtual void computeFlux
        (
            const label cellID,
            faceGaussField<vector>& F,
            const faceGaussField<scalar>& U
        ) = 0;

        virtual void computeFlux
        (
            const label cellID,
            faceGaussField<tensor>& F,
            const faceGaussField<vector>& U
        ) = 0;

        //- Reset if necessary
        virtual void reset() = 0;

    // Utility
        // - Construct orthonormal basis (n, t1, t2) from normal n
        //   (t1, t2 are tangential directions)
        void makeONB(const vector& n, vector& t1, vector& t2) const;

        // - Decompose vector U to normal and 2 tangential components
        void decomposeU
        (
            const vector& U,
            const vector& n,
            vector& Un,
            vector& Ut1,
            vector& Ut2
        ) const;

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
