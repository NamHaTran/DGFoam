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
    Foam::dgRoeFluxSolver

Description
    Roe approximate Riemann solver for DG methods.

SourceFiles
    dgRoeFluxSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgRoeFluxSolver_H
#define Foam_dgRoeFluxSolver_H

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
                    Class dgRoeFluxSolver Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgRoeFluxSolver
 * \brief Roe approximate Riemann solver for compressible Euler fluxes
 *
 * \details
 * This solver evaluates the Roe flux in a face-local orthonormal frame
 * whose first axis is the outward face normal. The left and right states are
 * projected into that frame, the one-dimensional Roe linearization is applied
 * to the full five-equation Euler system, and the resulting momentum flux is
 * rotated back to Cartesian coordinates.
 *
 * Supported transport equations:
 * - massTransport
 * - momentumTransport
 * - energyTransport
 *
 * Generic scalarTransport is intentionally rejected because the Roe solver
 * implemented here is an Euler-system flux, not a passive-scalar upwinder.
 */
class dgRoeFluxSolver
:
    public dgFluxSolver
{
private:

    /** Thermodynamics package used to recover interface gamma. */
    const dgThermoConservative& thermo_;

    /** Conservative density field. */
    const dgField<scalar>& rho_;

    /** Conservative total energy field. */
    const dgField<scalar>& E_;

    /** Velocity field used to build primitive states at the interface. */
    const dgField<vector>& U_;

    /** Pressure field used in the physical Euler fluxes. */
    const dgField<scalar>& p_;

    /** Emit the variable-gamma warning only once. */
    mutable bool warnedVariableGamma_;

    /** Assemble gamma at all Gauss points of the current cell. */
    void calcGamma
    (
        const label cellID,
        GaussField<scalar>& gamma
    ) const;

    /** Select the interface gamma used by the Roe linearization. */
    scalar selectInterfaceGamma
    (
        const scalar gammaL,
        const scalar gammaR
    ) const;

    /** Compute Roe mass, momentum, and energy fluxes at one face Gauss point. */
    void calcRoeFlux
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

    TypeName("Roe");

    dgRoeFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    virtual ~dgRoeFluxSolver() = default;

    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    virtual void reset() override {}

    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
