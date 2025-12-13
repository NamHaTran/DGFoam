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
    Foam::dgHLLEFluxSolver

Description
    HLLE (Harten–Lax–van Leer–Einfeldt) approximate Riemann solver.

    This flux uses only the bounding wave speeds S_L and S_R to construct a
    single intermediate state between the left (Minus/P-side) and right
    (Plus/N-side) states. HLLE does not resolve the contact or shear waves,
    and therefore produces a more diffusive but very robust and
    positivity-preserving flux.

    Key properties
    - Robust for strong shocks and large discontinuities.
    - Positivity-preserving for density and pressure under standard CFL
      conditions.
    - No use of a contact wave → more diffusive than HLLC.
    - No star-state reconstruction → simpler and cheaper to evaluate.
    - Suitable for any conserved quantity (rho, rhoE, rhoU, scalar).

    Procedure (local 1D Riemann problem along face normal n):
    1. Project left and right physical fluxes onto n.
    2. Compute left/right normal velocities and sound speeds.
    3. Estimate wave speeds S_L and S_R using either:
         * Davis (S_L = min(Un^- - a^-, Un^+ - a^+),
                  S_R = max(Un^- + a^-, Un^+ + a^+)), or
         * Roe–Einfeldt averaging for reduced numerical diffusion.
    4. Construct the HLLE flux:
          f_HLLE = f_L,                             if S_L >= 0
                 = f_R,                             if S_R <= 0
                 = ( S_R f_L - S_L f_R
                     + S_L S_R (U_R - U_L) ) / (S_R - S_L), otherwise
    5. For scalar equations the Cartesian flux is f_HLLE * n.
       For vector/momentum equations the Cartesian flux is f_HLLE as a vector
       obtained from the tensor flux projection; no multiplication by n.

    Notation
    - Minus  = owner cell (P-side)
    - Plus   = neighbour cell (N-side)
    - n      = unit normal pointing P → N
    - U_L/R  = conserved variables at Gauss point
    - f_L/R  = physical normal fluxes projected onto n

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
                           Class dgHLLEFluxSolver
\*---------------------------------------------------------------------------*/

class dgHLLEFluxSolver : public dgFluxSolver
{
public:

    //- Wave-speed estimate type
    enum speedEstimateType
    {
        seDavis,        //!< Side speeds: Un ± a
        seRoeEinfeldt   //!< Roe-averaged speeds
    };

private:

    //- Primitive fields (references)
        const dgField<scalar>& rho_;     //!< Density
        const dgField<vector>& U_;       //!< Velocity
        const dgField<scalar>& p_;       //!< Pressure
        const dgField<scalar>& a_;       //!< Speed of sound
        const dgField<scalar>& h_;       //!< Enthalpy
        const dgField<scalar>& gamma_;   //!< Cp/Cv

    //- Wavespeed selection
    speedEstimateType speedEst_;


    //- Compute S_L and S_R at Gauss point
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

    //- Runtime type info
    TypeName("HLLE");

    //- Constructor
    dgHLLEFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //- Destructor
    virtual ~dgHLLEFluxSolver() = default;

    //- Scalar U, vector F
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    //- Vector U, tensor F
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    //- HLLE requires no caching → empty
    virtual void reset() override {}

    //- Read settings
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
