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
    Foam::dgGeneralPDETerm

Description
    General DG PDE term for residual assembly. Handles operators of the form:

        Ri = - ∫ grad(phi_i) : F dV
             + ∫ phi_i (F_flux · n) dS

    where:
        - UType = scalar  and FType = vector
        - UType = vector  and FType = tensor

    F may be convective, diffusive, reacting, or any flux supplied by the user.

SourceFiles
    dgGeneralPDETerm.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeneralPDETerm_H
#define Foam_dgGeneralPDETerm_H

#include "word.H"
#include "List.H"
#include "label.H"
#include "pTraits.H"

#include "dgGeomMesh.H"
#include "dgBasisField.H"
#include "GaussField.H"
#include "dgFluxSolverManager.H"
#include "dgMath.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dgGeneralPDETerm Declaration
\*---------------------------------------------------------------------------*/

template<class UType, class FType>
class dgGeneralPDETerm
{
private:

    //- Name of this PDE term (used to fetch correct flux solver)
    word termName_;

    //- Equation type
    dgFluxSolver::equationType eqnType_;

    //- Cell index
    label cellID_;

    //- Reference to DG mesh
    dgGeomMesh& dgMesh_;

    //- Basis functions for this cell
    dgBasisField basis_;

    //- Solution field
    const GaussField<UType>& U_;

    //- Flux field (will be overwritten by flux solvers)
    GaussField<FType>& F_;

    //- Flux solver manager
    dgFluxSolverManager& fluxSolver_;

    //- Number of DoFs in this cell
    label nDof_;

    //- Residual per DoF
    List<UType> R_;


    //- Compile-time type check
    inline void checkTypes() const;

public:

    // Constructors

    //- Construct from mesh, term name, and fields
    dgGeneralPDETerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        const label cellID,
        dgGeomMesh& mesh,
        const GaussField<UType>& U,
        GaussField<FType>& F,
        dgFluxSolverManager& fluxSolver
    );

    //- Destructor
    ~dgGeneralPDETerm() = default;


    // Member Functions

    //- Return number of DoFs
    inline label nDof() const { return nDof_; }

    //- Return residual (const)
    inline const List<UType>& R() const { return R_; }

    //- Return residual (mutable)
    inline List<UType>& R() { return R_; }

    //- Compute numerical flux using the dgFluxSolverManager
    void updateFlux();

    //- Compute volume + surface contributions to residual
    void calcResidual();
};


} // End namespace Foam

#ifdef NoRepository
    #include "dgGeneralPDETermI.H"
#endif

#endif

// ************************************************************************* //
