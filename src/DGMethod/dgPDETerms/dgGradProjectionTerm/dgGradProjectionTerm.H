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
    Foam::dgGradProjectionTerm

Description
    DG auxiliary projection term for equations that project the gradient of
    a state field C into an auxiliary unknown.

        - volInt(C*grad(phi) dv) + surfInt(phi*flux(C) dS)

    This class stores/exposes the projected state field.  The residual is
    assembled by auxEqnResidual.  The face flux is not obtained from
    dgFluxSolverManager, and termName_ is not looked up in dgSchemes to select
    a flux solver. For auxiliary equations the residual assembler applies the
    default central face state:

        C* = 0.5*(C- + C+) on internal faces

    directly during surface integration. Boundary faces use the existing
    mid-state supplied by the boundary update path.

    Supported state/result pairs:
    - C scalar -> residual vector
    - C vector -> residual tensor with rows as vector components and columns
      as spatial directions

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGradProjectionTerm_H
#define Foam_dgGradProjectionTerm_H

#include "word.H"
#include "List.H"
#include "label.H"
#include "pTraits.H"
#include "products.H"

#include "dgBasisField.H"
#include "GaussField.H"
#include "dgFluxSolver.H"
#include "dgMath.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgGradProjectionTerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGradProjectionTerm
 * \brief Cell-local auxiliary gradient-projection term.
 *
 * \tparam UType Type of the projected state C. Supported: scalar or vector.
 *
 * Persistent GaussField inputs are referenced directly; temporaries are held
 * through tmp. This class does not assemble a residual by itself.  It exposes
 * the state field \f$C\f$ to auxEqnResidual, which evaluates
 * \f$-\int C\nabla\phi\,dV + \int \phi C^\ast n\,dS\f$.
 */
template<class UType>
class dgGradProjectionTerm
{
private:

    /** Name of this projection term; not used for dgSchemes flux lookup. */
    word termName_;

    /** Equation type tag. Defaults to auxiliaryEqn. */
    dgFluxSolver::equationType eqnType_;

    /** Owner cell index. */
    label cellID_;

    /** Basis functions associated with the owner cell. */
    const dgBasisField& basis_;

    /** Projected state C evaluated at cell and face Gauss points. */
    tmp<GaussField<UType>> U_;

    /** Compile-time validation for supported state/result pairs. */
    inline void checkTypes() const;

public:

    /** Type of the projected state C. */
    using state_type = UType;

    /** Type of C*grad(phi) and C*n. */
    using residual_type = typename outerProduct<UType, vector>::type;

    /**
     * \brief Construct from a persistent GaussField state C.
     *
     * The state is referenced directly; central internal-face states are
     * evaluated on read during calcResidual().
     */
    dgGradProjectionTerm
    (
        const word& termName,
        const dgBasisField& basis,
        const GaussField<UType>& U,
        const dgFluxSolver::equationType eqnType =
            dgFluxSolver::equationType::auxiliaryEqn
    );

    /**
     * \brief Construct from a temporary GaussField state C.
     *
     * The temporary is consumed when possible to avoid an extra copy.
     */
    dgGradProjectionTerm
    (
        const word& termName,
        const dgBasisField& basis,
        const tmp<GaussField<UType>>& tU,
        const dgFluxSolver::equationType eqnType =
            dgFluxSolver::equationType::auxiliaryEqn
    );

    ~dgGradProjectionTerm() = default;

    inline const GaussField<UType>& stateField() const { return U_(); }

    inline const GaussField<UType>* stateFieldPtr() const { return &U_(); }

    inline dgFluxSolver::equationType equationType() const { return eqnType_; }

    /**
     * \brief Compatibility hook for old callers.
     *
     * Central internal-face states are evaluated directly in auxEqnResidual,
     * so this function does not mutate faceField().midValues().
     */
    void updateFlux();
};

} // End namespace Foam

#ifdef NoRepository
    #include "dgGradProjectionTermI.H"
#endif

#endif

// ************************************************************************* //
