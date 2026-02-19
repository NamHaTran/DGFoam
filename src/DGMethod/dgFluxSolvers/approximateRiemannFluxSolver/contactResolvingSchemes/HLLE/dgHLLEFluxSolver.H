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
    Foam::dgHLLEFluxSolver

Description
    HLLE approximate Riemann solver for DG methods.

SourceFiles
    dgHLLEFluxSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgHLLEFluxSolver_H
#define Foam_dgHLLEFluxSolver_H

#include "dgFluxSolver.H"
#include "dgField.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

template<class Type> class faceGaussField;

/*---------------------------------------------------------------------------*\
                     Class dgHLLEFluxSolver Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgHLLEFluxSolver
 * \brief HLLE (Harten–Lax–van Leer–Einfeldt) numerical flux
 *
 * \details
 * Implements the HLLE approximate Riemann solver for Discontinuous
 * Galerkin methods.
 *
 * HLLE resolves the local one-dimensional Riemann problem using only
 * the bounding wave speeds \f$S_L\f$ and \f$S_R\f$, without restoring
 * the contact or shear waves.
 *
 * The numerical flux is defined as:
 * \f[
 * \boldsymbol{f}_{\mathrm{HLLE}} =
 * \begin{cases}
 * \boldsymbol{f}_L,
 *   & S_L \ge 0, \\
 * 
 * \boldsymbol{f}_R,
 *   & S_R \le 0, \\
 * 
 * \frac{
 *   S_R \boldsymbol{f}_L
 * - S_L \boldsymbol{f}_R
 * + S_L S_R (\boldsymbol{U}_R - \boldsymbol{U}_L)
 * }{S_R - S_L},
 *   & \text{otherwise}
 * \end{cases}
 * \f]
 *
 * where:
 * - \f$\boldsymbol{U}_L, \boldsymbol{U}_R\f$ are conserved variables
 * - \f$\boldsymbol{f}_L, \boldsymbol{f}_R\f$ are physical fluxes projected
 *   onto the face normal
 *
 * Wave speeds are estimated using either:
 * - Davis estimate: \f$U_n \pm a\f$
 * - Roe–Einfeldt averaged estimate
 *
 * Properties:
 * - Robust for strong shocks
 * - Positivity-preserving under standard CFL conditions
 * - More diffusive than HLLC
 * - No star-state reconstruction
 *
 * Cartesian flux handling:
 * - Scalar transport: returns \f$f_n \boldsymbol{n}\f$
 * - Vector/momentum transport: returns vector flux directly
 */
class dgHLLEFluxSolver
:
    public dgFluxSolver
{
public:

    /** Wave-speed estimation method */
    enum speedEstimateType
    {
        seDavis,        //!< Davis side-speed estimate
        seRoeEinfeldt   //!< Roe–Einfeldt averaged estimate
    };

private:

    /** Density field (read-only, non-owning) */
    const dgField<scalar>& rho_;

    /** Velocity field (read-only, non-owning) */
    const dgField<vector>& U_;

    /** Pressure field (read-only, non-owning) */
    const dgField<scalar>& p_;

    /** Speed of sound field (read-only, non-owning) */
    const dgField<scalar>& a_;

    /** Enthalpy field (read-only, non-owning) */
    const dgField<scalar>& h_;

    /** Heat-capacity ratio Cp/Cv (read-only, non-owning) */
    const dgField<scalar>& gamma_;

    /** Selected wave-speed estimation method */
    speedEstimateType speedEst_;

    /**
     * \brief Compute HLLE wave speeds
     *
     * Computes the left and right bounding wave speeds
     * \f$S_L\f$ and \f$S_R\f$ at a face Gauss point.
     */
    void calcWaveSpeed
    (
        scalar& SL,
        scalar& SR,
        const scalar rhoL,
        const scalar rhoR,
        const vector& ULv,
        const vector& URv,
        const scalar pL,
        const scalar pR,
        const scalar aL,
        const scalar aR,
        const scalar gammaL,
        const scalar gammaR,
        const scalar hL,
        const scalar hR,
        const vector& n
    ) const;

public:

    //- Runtime type information
    TypeName("HLLE");

    /**
     * \brief Construct HLLE flux solver
     *
     * Reads solver settings and initializes references
     * to required flow fields.
     */
    dgHLLEFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** Destructor */
    virtual ~dgHLLEFluxSolver() = default;

    /**
     * \brief Compute numerical flux for scalar conserved variable
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    /**
     * \brief Compute numerical flux for vector conserved variable
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    /** HLLE has no internal cache */
    virtual void reset() override {}

    /**
     * \brief Read solver configuration
     */
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
