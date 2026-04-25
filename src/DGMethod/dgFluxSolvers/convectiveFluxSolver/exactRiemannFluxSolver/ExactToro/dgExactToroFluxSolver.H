/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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
    Foam::dgExactToroFluxSolver

Description
    Exact Euler Riemann solver based on Toro (2009).

SourceFiles
    dgExactToroFluxSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgExactToroFluxSolver_H
#define Foam_dgExactToroFluxSolver_H

#include "dgFluxSolver.H"
#include "dgField.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{
    class dgThermoConservative;
    template<class Type> class faceGaussField;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dgExactToroFluxSolver Declaration
\*---------------------------------------------------------------------------*/

/**
 * \brief Exact numerical flux for the compressible Euler equations using
 * Toro's exact 1-D Riemann solver.
 *
 * This solver evaluates the interfacial flux by solving the local
 * left/right Riemann problem in the face-normal direction. The full
 * three-dimensional velocity is first projected onto a local orthonormal
 * basis \f$(\mathbf{n}, \mathbf{t}_1, \mathbf{t}_2)\f$, so that the exact
 * solution procedure is applied only to the normal velocity component while
 * the tangential components are transported through the sampled star state.
 *
 * The implementation follows the standard Toro workflow:
 * - build primitive left/right states at each face Gauss point,
 * - estimate the star-region pressure with Toro's pressure guess,
 * - refine \f$p^*\f$ by Newton iteration on the exact pressure function,
 * - recover \f$u^*\f$ and sample the self-similar solution at
 *   \f$\xi = x/t = 0\f$,
 * - reconstruct mass, momentum, and energy fluxes from the sampled state.
 *
 * \note The exact solver assumes a calorically perfect gas with constant
 * \f$\gamma\f$. If the left and right states provide different values, this
 * implementation emits a warning once and uses the arithmetic mean at the
 * interface.
 */
class dgExactToroFluxSolver
:
    public dgFluxSolver
{
private:

    /// Thermodynamic package used to recover temperature-dependent properties.
    const dgThermoConservative& thermo_;

    /// Conservative density field.
    const dgField<scalar>& rho_;

    /// Conservative total-energy field.
    const dgField<scalar>& E_;

    /// Velocity field used to form primitive states at the interface.
    const dgField<vector>& U_;

    /// Thermodynamic pressure field.
    const dgField<scalar>& p_;

    /// Ensures that the variable-gamma warning is emitted only once.
    mutable bool warnedVariableGamma_;

    /// Ensures that the vacuum-fallback warning is emitted only once.
    mutable bool warnedVacuumFallback_;

    /// Ensures that the non-convergence fallback warning is emitted only once.
    mutable bool warnedNewtonFallback_;

    /**
     * \brief Compute the ratio of specific heats at all Gauss points of a cell.
     * \param cellID Cell for which the Gauss field is assembled.
     * \param gamma Output Gauss field of \f$\gamma = C_p/C_v\f$.
     */
    void calcGamma
    (
        const label cellID,
        GaussField<scalar>& gamma
    ) const;

    /**
     * \brief Select the interface value of \f$\gamma\f$ used by the exact solver.
     * \param gammaL Left-state ratio of specific heats.
     * \param gammaR Right-state ratio of specific heats.
     * \return Interface value of \f$\gamma\f$ used in the Riemann solve.
     */
    scalar selectInterfaceGamma
    (
        const scalar gammaL,
        const scalar gammaR
    ) const;

    /**
     * \brief Solve the local exact Riemann problem and return Euler fluxes.
     *
     * The input states are interpreted as left/right primitive data at a
     * single face Gauss point. The routine rotates the velocity into the
     * face-local frame, solves the 1-D normal Riemann problem, samples the
     * state at the interface, and reconstructs the physical flux vector.
     *
     * \param rhoL Left density.
     * \param rhoR Right density.
     * \param ULv Left velocity vector.
     * \param URv Right velocity vector.
     * \param EL Left total energy density.
     * \param ER Right total energy density.
     * \param pL Left pressure.
     * \param pR Right pressure.
     * \param gammaL Left ratio of specific heats.
     * \param gammaR Right ratio of specific heats.
     * \param n Unit face normal pointing from minus to plus side.
     * \param massFlux Output mass flux \f$\rho u_n\f$.
     * \param momentumFlux Output momentum flux tensor row projected onto
     * \f$\mathbf{n}\f$.
     * \param energyFlux Output total-energy flux.
     */
    void calcExactFlux
    (
        const scalar rhoL,
        const scalar rhoR,
        const vector& ULv,
        const vector& URv,
        const scalar EL,
        const scalar ER,
        const scalar pL,
        const scalar pR,
        const scalar gammaL,
        const scalar gammaR,
        const vector& n,
        scalar& massFlux,
        vector& momentumFlux,
        scalar& energyFlux
    ) const;

public:

    TypeName("ExactToro");

    /// Construct from run-time dictionary and mesh context.
    dgExactToroFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    virtual ~dgExactToroFluxSolver() = default;

    /**
     * \brief Compute scalar-valued fluxes.
     *
     * This overload is used for mass and total-energy transport equations.
     * The exact Riemann problem is solved independently at each face Gauss
     * point and the resulting scalar flux is projected along the face normal.
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    /**
     * \brief Compute vector-valued fluxes for the momentum equation.
     *
     * The returned tensor-valued face flux is the exact momentum flux sampled
     * from the interfacial Riemann solution.
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    /// No cached state is stored between calls.
    virtual void reset() override {}

    /// Read optional solver controls from the dictionary.
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
