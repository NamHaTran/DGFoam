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

Description
    Central numerical flux for DG methods.

SourceFiles
    dgCentralFluxSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCentralFluxSolver_H
#define Foam_dgCentralFluxSolver_H

#include "dgFluxSolver.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{
    template<class Type> class faceGaussField;
}

/*---------------------------------------------------------------------------*\
                  Class dgCentralFluxSolver Declaration
\*---------------------------------------------------------------------------*/

namespace Foam
{

/**
 * \class Foam::dgCentralFluxSolver
 * \brief Central (non-dissipative) numerical flux
 *
 * \details
 * Implements the classical central flux used in DG formulations.
 *
 * The numerical flux is defined as the arithmetic average of the
 * physical flux evaluated on the minus and plus sides:
 *
 * \f[
 *   \boldsymbol{h}
 *   =
 *   \frac{1}{2}
 *   \left(
 *     \boldsymbol{F}^{-} + \boldsymbol{F}^{+}
 *   \right)
 * \f]
 *
 * where:
 * - \f$\boldsymbol{F}^{-}\f$ is the flux from the interior (minus-side) state
 * - \f$\boldsymbol{F}^{+}\f$ is the flux from the exterior (plus-side) state
 *
 * No upwinding or numerical dissipation is applied. This flux is therefore
 * non-dissipative and non-directional.
 *
 * Typical use cases:
 * - Diffusive fluxes
 * - LDG auxiliary equations
 * - Elliptic operators with external stabilisation
 *
 * This flux is generally unsuitable for convection-dominated hyperbolic
 * problems without additional stabilisation.
 */
class dgCentralFluxSolver
:
    public dgFluxSolver
{
public:

    /** Runtime type information */
    TypeName("central");

    // Constructors --------------------------------------------------------

    /**
     * \brief Construct central flux solver
     *
     * \param name Solver name
     * \param dict Flux-solver dictionary
     * \param mesh DG geometric mesh
     */
    dgCentralFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    )
    :
        dgFluxSolver(name, dict, mesh)
    {}

    //- Destructor
    virtual ~dgCentralFluxSolver() = default;

    // Flux evaluation -----------------------------------------------------

    /**
     * \brief Compute central flux for scalar field
     *
     * Computes the arithmetic average of vector-valued fluxes
     * corresponding to a scalar conserved variable.
     *
     * \param cellID DG cell index
     * \param F      Flux field (output)
     * \param U      Scalar solution field at face Gauss points
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    /**
     * \brief Compute central flux for vector field
     *
     * Computes the arithmetic average of tensor-valued fluxes
     * corresponding to vector-valued conserved variables.
     *
     * \param cellID DG cell index
     * \param F      Flux field (output)
     * \param U      Vector solution field at face Gauss points
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    // Optional hooks ------------------------------------------------------

    /** No internal state to reset */
    virtual void reset() override {}

    /** No solver-specific parameters to read */
    virtual void read(const dictionary& dict) override {}
};

} // End namespace Foam

#endif
// ************************************************************************* //
