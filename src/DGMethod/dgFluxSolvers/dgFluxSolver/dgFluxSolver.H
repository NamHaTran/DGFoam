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
    Base class for DG numerical flux solvers.

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

/**
 * \enum equationType
 * \class Foam::dgFluxSolver
 * \brief Base class for DG numerical flux solvers
 *
 * \details
 * Defines the runtime-selectable interface for numerical-flux solvers that
 * operate on faceGaussField data. Implementations compute a numerical flux
 * from minus/plus states at face Gauss points.
 *
 * Plus/minus convention:
 * - minus: owner cell (P-side)
 * - plus : neighbour cell (N-side)
 *
 * The equationType is used by derived solvers to interpret the meaning of
 * state components and to select wave-speed estimates where required.
 */
class dgFluxSolver
{
public:
    // Enumeration for identifying the PDE being solved
    enum class equationType
    {
        /** Mass transport equation */
        massTransport,

        /** Momentum transport equation */
        momentumTransport,

        /** Energy transport equation */
        energyTransport,

        /** Scalar transport equation */
        scalarTransport
    };

protected:

    // Attributes

        /** Solver instance name (for logging/debug) */
        word name_;

        /** Configuration dictionary for this solver */
        dictionary dict_;

        /** Reference to DG geometric mesh */
        const dgGeomMesh& mesh_;

        /** Equation type used for interpreting left/right states */
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

        /**
         * \brief Construct from name and dictionary
         *
         * \param name Solver instance name
         * \param dict Solver dictionary
         * \param mesh DG geometric mesh
         */
        dgFluxSolver(const word& name, const dictionary& dict, const dgGeomMesh& mesh);

        //- Disallow default construction
        dgFluxSolver() = delete;

        //- Destructor
        virtual ~dgFluxSolver() = default;

    // Factory

        /**
         * \brief Select and construct from dictionary
         *
         * \param name           Solver instance name
         * \param fluxSolverType Selected flux-solver type name
         * \param dict           Solver dictionary
         * \param mesh           DG geometric mesh
         * \return autoPtr to constructed dgFluxSolver
         */
        static autoPtr<dgFluxSolver> New
        (
            const word& name,
            const word& fluxSolverType,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

    // Access

        /** Return solver name */
        inline const word& name() const { return name_; }

        /** Return solver dictionary */
        inline const dictionary& dict() const { return dict_; }

    // Core API

        /**
         * \brief Compute numerical flux for scalar conserved variable
         *
         * \param cellID DG cell index
         * \param F      Flux field (output)
         * \param U      Scalar solution at face Gauss points
         */
        virtual void computeFlux
        (
            const label cellID,
            faceGaussField<vector>& F,
            const faceGaussField<scalar>& U
        ) = 0;

        /**
         * \brief Compute numerical flux for vector conserved variable
         *
         * \param cellID DG cell index
         * \param F      Flux field (output)
         * \param U      Vector solution at face Gauss points
         */
        virtual void computeFlux
        (
            const label cellID,
            faceGaussField<tensor>& F,
            const faceGaussField<vector>& U
        ) = 0;

        /** Reset solver internal cache/state if required */
        virtual void reset() = 0;

    // Utility

        /**
         * \brief Construct an orthonormal basis from a normal vector
         *
         * Builds tangential directions (t1, t2) such that (n, t1, t2) is an
         * orthonormal triad.
         *
         * \param n  Unit normal
         * \param t1 First tangential direction (output)
         * \param t2 Second tangential direction (output)
         */
        void makeONB(const vector& n, vector& t1, vector& t2) const;

        /**
         * \brief Decompose a vector into normal and tangential components
         *
         * \param U   Input vector
         * \param n   Unit normal
         * \param Un  Normal component (output)
         * \param Ut1 Tangential component along t1 (output)
         * \param Ut2 Tangential component along t2 (output)
         */
        void decomposeU
        (
            const vector& U,
            const vector& n,
            vector& Un,
            vector& Ut1,
            vector& Ut2
        ) const;

        /**
         * \brief Read solver settings from dictionary
         *
         * \param dict Solver dictionary
         */
        virtual void read(const dictionary& dict) = 0;

        /** Set equation type */
        inline void setEquationType(equationType t) { eqnType_ = t; }

        /** Return equation type */
        inline equationType getEquationType() const { return eqnType_; }
};

    /**
     * \brief Stream output for dgFluxSolver::equationType
     *
     * \param os   Output stream
     * \param type Equation type
     * \return Output stream
     */
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
