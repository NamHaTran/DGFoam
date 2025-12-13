/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgHLLCFluxSolver

Description
    HLLC (Harten–Lax–van Leer–Contact) approximate Riemann solver.

    This flux restores the missing contact and shear waves that are removed in
    the HLLE formulation. HLLC introduces a middle wave S* and two
    reconstructed intermediate states (U*_L and U*_R), allowing accurate
    resolution of contact discontinuities and tangential shear layers.

    Compared to HLLE:
    - HLLC is significantly less diffusive.
    - Captures contact discontinuities sharply.
    - Preserves tangential velocity continuity (important for momentum flux).
    - Requires computing a star region and wave strengths C_L, C_R.

    Procedure (local 1D Riemann problem along face normal n):
    1. Project left and right physical fluxes onto n.
    2. Compute left/right normal velocities and sound speeds.
    3. Estimate wave speeds S_L and S_R using either:
         * Davis side speeds, or
         * Roe–Einfeldt averages (U_roe, H_roe, a_roe).
    4. Compute the contact wave speed S* based on pressure and velocity jumps.
    5. Compute the scaling coefficients:
         C_L = rho_L (S_L - Un_L) / (S_L - S*)
         C_R = rho_R (S_R - Un_R) / (S_R - S*)
    6. Construct the intermediate star states U*_L and U*_R:
         - For mass:     U* = C
         - For momentum: U* = C ( S* n + U_t ), preserving tangential velocity
         - For energy:   U* = C ( e + corrective term )
         - For scalars:  U* = C q
    7. Construct the flux in each region:
         f = f_L,                               if S_L >= 0
         f = f*_L,                              if S_L <= 0 <= S*
         f = f*_R,                              if S* <= 0 <= S_R
         f = f_R,                               if S_R <= 0
       where f*_L and f*_R are star-region fluxes.
    8. Cartesian flux:
         * scalar equations → f * n
         * momentum equations → f as a vector (no multiplication by n)

    Notation
    - Minus  = owner cell (P-side)
    - Plus   = neighbour cell (N-side)
    - n      = unit normal pointing P → N
    - U_L/R  = conserved variables at Gauss point
    - f_L/R  = physical normal fluxes projected onto n

    Advantages
    - Accurate contact resolution.
    - Good shear preservation.
    - Less diffusive than HLLE while retaining robustness.

    References
    - Toro (2009), “Riemann Solvers and Numerical Methods for Fluid Dynamics”
    - Harten, Lax, van Leer (1983)
    - Einfeldt (1991)

SourceFiles
    dgHLLCFluxSolver.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_dgHLLCFluxSolver_H
#define Foam_dgHLLCFluxSolver_H

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
                        Class dgHLLCFluxSolver Declaration
\*---------------------------------------------------------------------------*/

class dgHLLCFluxSolver : public dgFluxSolver
{
public:

    //- Wave-speed estimate types
    enum speedEstimateType
    {
        seDavis,        // S_L = min(Un^- - a^-, Un^+ - a^+), S_R = max(Un^- + a^-, Un^+ + a^+)
        seRoeEinfeldt   // Roe averages (U_roe, H_roe, a_roe)
    };

private:
    // Attributes

        //- Reference to global fields
            // Density
            const dgField<scalar>& rho_;
            // Velocity
            const dgField<vector>& U_;
            // Pressure
            const dgField<scalar>& p_;
            // Speed of sound
            const dgField<scalar>& a_;
            // Enthalpy
            const dgField<scalar>& h_;
            // Ratio Cp by Cv
            const dgField<scalar>& gamma_;

        //- Wave-speed estimate type
        speedEstimateType speedEst_;

        //- List intermediate states
        List<List<scalar>> SL_list_;
        List<List<scalar>> SR_list_;
        List<List<scalar>> SStar_list_;
        List<List<scalar>> CL_list_;
        List<List<scalar>> CR_list_;

        //- List of computed intermediate states flag
        List<bool> isStateComputed_;

public:

    //- Runtime type name
    TypeName("HLLC");

    //- Construct from name, dictionary and mesh
    dgHLLCFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //- Destructor
    virtual ~dgHLLCFluxSolver() = default;

    //- Calculate intermediate states
    void calcIntermediateState
    (
        const label cellID,
        const label localFaceID,
        const label localGaussID,
        const vector& n,

        // Output
        scalar& SL,
        scalar& SR,
        scalar& SStar,
        scalar& CL,
        scalar& CR
    );

    //- Compute numerical flux from Plus/Minus physical flux vectors
    //  Input:
    //    - cellID  : cellID of cell where flux is calculated
    //    - F       : faceGaussField of field F
    //    - U       : faceGaussField of U
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

    // Override reset() method
    virtual void reset() override;

    //- Reading method
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
