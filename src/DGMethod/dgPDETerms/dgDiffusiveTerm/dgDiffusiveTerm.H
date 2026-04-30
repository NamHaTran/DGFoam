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
    Foam::dgDiffusiveTerm

Description
    Diffusive DG PDE term used to compute and expose fluxes.

SourceFiles
    dgDiffusiveTerm.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgDiffusiveTerm_H
#define Foam_dgDiffusiveTerm_H

#include "word.H"
#include "List.H"
#include "label.H"
#include "pTraits.H"

#include "dgBasisField.H"
#include "GaussField.H"
#include "dgFluxSolverManager.H"
#include "dgMath.H"
#include "dgPDETerm.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dgDiffusiveTerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgDiffusiveTerm
 * \brief Diffusive DG PDE term for flux assembly.
 *
 * \tparam UType Type of solution variable (scalar or vector)
 * \tparam FType Type of diffusive flux (vector or tensor)
 *
 * The constructor mirrors dgConvectiveTerm. The numerical surface flux is
 * selected through dgFluxSolverManager using the PDE term name, and the
 * resulting flux field is integrated by dgResidual.
 */
template<class UType, class FType>
class dgDiffusiveTerm
:
    public dgPDETerm<UType, FType>
{
private:

    /** Name of PDE term, kept for symmetry with dgConvectiveTerm. */
    word termName_;

    /** Equation type tag, kept for symmetry with dgConvectiveTerm. */
    dgFluxSolver::equationType eqnType_;

    /** Flux type required by this PDE term. */
    dgFluxSolver::fluxType fluxType_;

    /** Owner cell index. */
    label cellID_;

    /** Basis functions associated with this cell. */
    const dgBasisField& basis_;

    /** Solution field evaluated at Gauss points. */
    const GaussField<UType>& U_;

    /** Diffusive flux field evaluated at Gauss points. */
    tmp<GaussField<FType>> F_;

    /** Flux solver manager. */
    dgFluxSolverManager& fluxSolver_;

    /** Flux solver selected for this term. */
    dgFluxSolver& termFluxSolver_;

    /** Compile-time validation of UType/FType combination. */
    inline void checkTypes() const;

public:

    using solution_type = UType;
    using flux_type = FType;

    dgDiffusiveTerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        dgFluxSolver::fluxType fluxType,
        const dgBasisField& basis,
        const GaussField<UType>& U,
        GaussField<FType>& F,
        dgFluxSolverManager& fluxSolver
    );

    dgDiffusiveTerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        dgFluxSolver::fluxType fluxType,
        const dgBasisField& basis,
        const GaussField<UType>& U,
        const tmp<GaussField<FType>>& tF,
        dgFluxSolverManager& fluxSolver
    );

    ~dgDiffusiveTerm() = default;

    inline bool isDiffusion() const override { return true; }

    inline const GaussField<UType>& stateField() const { return U_; }

    inline const GaussField<FType>& fluxField() const { return F_(); }

    inline const GaussField<FType>* fluxFieldPtr() const override
    {
        return &F_();
    }

    inline const word& termName() const { return termName_; }

    inline dgFluxSolver::equationType equationType() const { return eqnType_; }

    inline dgFluxSolver::fluxType fluxType() const { return fluxType_; }

    /**
     * \brief Compute numerical surface flux.
     *
     * Selects the appropriate dgFluxSolver via dgFluxSolverManager
     * and updates the face Gauss flux field.
     */
    void updateFlux();

};

} // End namespace Foam

#ifdef NoRepository
    #include "dgDiffusiveTermI.H"
#endif

#endif

// ************************************************************************* //
