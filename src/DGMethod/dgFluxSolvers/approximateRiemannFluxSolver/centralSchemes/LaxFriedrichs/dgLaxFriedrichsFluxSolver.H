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
    Foam::dgLaxFriedrichsFluxSolver

Description
    Local Lax–Friedrichs (Rusanov) numerical flux for DG methods.

SourceFiles
    dgLaxFriedrichsFluxSolver.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_dgLaxFriedrichsFluxSolver_H
#define Foam_dgLaxFriedrichsFluxSolver_H

#include "dgFluxSolver.H"
#include "dgField.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{
    template<class Type> class faceGaussField;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dgLaxFriedrichsFluxSolver Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgLaxFriedrichsFluxSolver
 * \brief Local Lax–Friedrichs (Rusanov) numerical flux
 *
 * \details
 * Implements the local Lax–Friedrichs (Rusanov) flux for DG formulations.
 *
 * The numerical flux is defined as:
 * \f[
 *   \boldsymbol{h}
 *   =
 *   \frac{1}{2}
 *   \left(
 *     \boldsymbol{F}(\boldsymbol{U}^{-})
 *   + \boldsymbol{F}(\boldsymbol{U}^{+})
 *   \right)\cdot\boldsymbol{n}
 *   -
 *   \frac{1}{2} C
 *   \left(
 *     \boldsymbol{U}^{+} - \boldsymbol{U}^{-}
 *   \right)
 * \f]
 *
 * where:
 * - \f$\boldsymbol{U}^{-}\f$ is the interior (minus-side) state
 * - \f$\boldsymbol{U}^{+}\f$ is the exterior (plus-side) state
 * - \f$\boldsymbol{F}\f$ is the physical flux
 * - \f$\boldsymbol{n}\f$ is the outward unit normal
 *
 * The dissipation coefficient \f$C\f$ is evaluated locally from flow
 * velocity and acoustic speed. Optional Mach-based scaling can be enabled
 * via dictionary settings.
 *
 * This flux is robust but introduces stronger numerical dissipation than
 * approximate Riemann solvers such as HLLC or Roe.
 */
class dgLaxFriedrichsFluxSolver
:
    public dgFluxSolver
{
private:

    /** Pointer to velocity field (read-only, non-owning) */
    const dgField<vector>& U_;

    /** Pointer to speed-of-sound field (read-only, non-owning) */
    const dgField<scalar>& a_;

    /** Enable Mach-number-based dissipation scaling */
    bool scaleByMach_;

    /**
     * \brief Compute local dissipation coefficient
     *
     * Computes the Lax–Friedrichs dissipation coefficient using left/right
     * states and face normal direction.
     *
     * \param ULv Velocity on minus side
     * \param URv Velocity on plus side
     * \param aL  Speed of sound on minus side
     * \param aR  Speed of sound on plus side
     * \param n   Face outward unit normal
     * \return Dissipation coefficient
     */
    scalar calcDissipationCoeff
    (
        const vector& ULv,
        const vector& URv,
        const scalar aL,
        const scalar aR,
        const vector& n
    ) const;

public:

    //- Runtime type information
    TypeName("LaxFriedrichs");

    /**
     * \brief Construct from name, dictionary, and DG mesh
     *
     * Initializes field references and reads solver options
     * (e.g. Mach-based scaling).
     *
     * \param name Solver name
     * \param dict Flux-solver dictionary
     * \param mesh DG geometric mesh
     */
    dgLaxFriedrichsFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //- Destructor
    virtual ~dgLaxFriedrichsFluxSolver() = default;

    /**
     * \brief Compute numerical flux for scalar conserved variable
     *
     * Produces a vector-valued numerical flux.
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
    ) override;

    /**
     * \brief Compute numerical flux for vector conserved variable
     *
     * Produces a tensor-valued numerical flux.
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
    ) override;

    /** Reset internal state (no-op) */
    virtual void reset() override {}

    /**
     * \brief Read solver settings from dictionary
     *
     * Updates dissipation-scaling options.
     *
     * \param dict Flux-solver dictionary
     */
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif
// ************************************************************************* //
