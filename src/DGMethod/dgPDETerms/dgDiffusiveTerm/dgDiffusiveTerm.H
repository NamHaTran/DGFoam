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
    Diffusive DG PDE term used for residual assembly with a fixed central flux.

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
#include "dgFluxSolver.H"
#include "dgMath.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dgDiffusiveTerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgDiffusiveTerm
 * \brief Diffusive DG PDE term for residual assembly.
 *
 * \tparam UType Type of solution variable (scalar or vector)
 * \tparam FType Type of diffusive flux (vector or tensor)
 *
 * The constructor mirrors dgConvectiveTerm except that no
 * dgFluxSolverManager is needed. The numerical surface flux is always the
 * central flux:
 *
 *     F* . n = 0.5*(F- . n + F+ . n)
 *
 * Boundary faces use F+ = F-. This is the intended path for fluxes built from
 * gradient fields where the boundary exterior gradient is prescribed equal to
 * the interior gradient.
 */
template<class UType, class FType>
class dgDiffusiveTerm
{
private:

    /** Name of PDE term, kept for symmetry with dgConvectiveTerm. */
    word termName_;

    /** Equation type tag, kept for symmetry with dgConvectiveTerm. */
    dgFluxSolver::equationType eqnType_;

    /** Owner cell index. */
    label cellID_;

    /** Basis functions associated with this cell. */
    const dgBasisField& basis_;

    /** Solution field evaluated at Gauss points. */
    const GaussField<UType>& U_;

    /** Diffusive flux field evaluated at Gauss points. */
    tmp<GaussField<FType>> F_;

    /** Number of degrees of freedom in this cell. */
    label nDof_;

    /** Residual values per degree of freedom. */
    List<UType> R_;

    /** Compile-time validation of UType/FType combination. */
    inline void checkTypes() const;

    /** Apply the fixed central flux to scalar-equation vector fluxes. */
    static void updateCentralFlux(faceGaussField<vector>& F);

    /** Apply the fixed central flux to vector-equation tensor fluxes. */
    static void updateCentralFlux(faceGaussField<tensor>& F);

public:

    using solution_type = UType;
    using flux_type = FType;

    dgDiffusiveTerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        const dgBasisField& basis,
        const GaussField<UType>& U,
        GaussField<FType>& F
    );

    dgDiffusiveTerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        const dgBasisField& basis,
        const GaussField<UType>& U,
        const tmp<GaussField<FType>>& tF
    );

    ~dgDiffusiveTerm() = default;

    inline label nDof() const { return nDof_; }

    inline const List<UType>& R() const { return R_; }

    inline List<UType>& R() { return R_; }

    inline const GaussField<UType>& stateField() const { return U_; }

    inline const GaussField<FType>& fluxField() const { return F_(); }

    inline const word& termName() const { return termName_; }

    inline dgFluxSolver::equationType equationType() const { return eqnType_; }

    /**
     * \brief Compute central numerical surface flux.
     *
     * Physical boundary faces first copy plus-side flux values from
     * minus-side values. Internal and processor faces use their existing
     * plus/minus traces in the same central formula.
     */
    void updateFlux();

    /**
     * \brief Assemble cell residual with the same sign convention as
     * dgConvectiveTerm.
     */
    void calcResidual();
};

} // End namespace Foam

#ifdef NoRepository
    #include "dgDiffusiveTermI.H"
#endif

#endif

// ************************************************************************* //
