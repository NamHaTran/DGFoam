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
    Foam::dgRoeFluxSolver

Description
    HLLC (Harten–Lax–van Leer–Contact) numerical flux for a single
    conserved equation working on vector-valued physical fluxes.

    Conventions
    - Plus  = neighbour cell (N-side)
    - Minus = owner cell (P-side)
    - 'n' is the unit normal pointing from P to N.
    - 'gaussID' is the local Gauss point index on the face.
    - UR, UL are the conserved values on Plus (N) and Minus (P) sides.

    Inputs available at the face Gauss point
    - From fieldsContext (face Gauss fields):
      * rho^± : density (scalar),   U^± : velocity (vector),   p^± : pressure (scalar)
    - From thermo (fluid properties handles):
      * T^± = eqnOfState().T(rho^±, p^±)
      * a^± = thermo().a(T^±)                // speed of sound
      * gamma = thermo().gamma()             // ratio of specific heats (if needed)

    HLLC procedure (local 1D Riemann along n)
    1. Local projection to normal direction:
       - f_R = FPlus  · n        // right  (Plus, N-side)  physical flux dotted with n
       - f_L = FMinus · n        // left   (Minus, P-side)
       - U_n^± = U^± · n         // normal velocity on each side

    2. Acoustic speeds on each side using thermo:
       - T^± = eqnOfState().T(rho^±, p^±)
       - a^± = thermo().a(T^±)

    3. Wave-speed bounds (Davis or Roe-Einfeldt):
       - S_L = min( U_n^- - a^- , Roe-estimate )
       - S_R = max( U_n^+ + a^+ , Roe-estimate )

    4. Contact wave speed (middle wave):
       - S* = (p^+ - p^- + rho^-·U_n^-·(S_L - U_n^-) - rho^+·U_n^+·(S_R - U_n^+)) /
              (rho^-·(S_L - U_n^-) - rho^+·(S_R - U_n^+))

    5. Intermediate states (U*^±), depend on equationType:
       - C^± = rho^±·(S_{±} - U_n^±) / (S_{±} - S*)
       - U*_L, U*_R computed based on equation:
         e.g. momentum → C^±·S*, energy → C^±·(e + φ), scalar → C^±·q

    6. Intermediate fluxes (F*^±):
       - F*_L = FMinus + S_L·(U*_L - UL)
       - F*_R = FPlus  + S_R·(U*_R - UR)

    7. HLLC normal numerical flux f_n:
       - if S_L ≥ 0       : f_n = f_L
       - else if S* ≥ 0   : f_n = F*_L · n
       - else if S_R > 0  : f_n = F*_R · n
       - else             : f_n = f_R

    8. Rotate back to Cartesian:
       - flux = f_n * n      // tangential parts vanish in 1D Riemann solve

    Notes
    - Ensure 'n' is a unit vector and oriented from P to N.
    - HLLC restores contact and shear waves, and is less diffusive than HLLE.
    - Conserved variable (UR/UL) must match the PDE being solved:
      (e.g. rho, rhoU, rhoE, or scalar transported quantity).

    References
    - Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics", 3rd ed., 2009.
    - Harten, Lax, van Leer (1983), Einfeldt (1991)
    - Hesthaven & Warburton, "Nodal Discontinuous Galerkin Methods", 2008.

SourceFiles
    dgRoeFluxSolver.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_dgRoeFluxSolver_H
#define Foam_dgRoeFluxSolver_H

#include "dgFluxSolver.H"
#include "vector.H"

namespace Foam
{
    template<class Type> class faceGaussField;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dgRoeFluxSolver Declaration
\*---------------------------------------------------------------------------*/

class dgRoeFluxSolver : public dgFluxSolver
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

        //- Wave-speed estimate type
        speedEstimateType speedEst_;

public:

    //- Runtime type name
    TypeName("HLLE");

    //- Construct from name, dictionary and thermo
    dgRoeFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgThermo& thermo
    );

    //- Destructor
    virtual ~dgRoeFluxSolver() = default;

    //- Compute numerical flux from Plus/Minus physical flux vectors
    //  Input:
    //    - gaussID : local Gauss point ID on face
    //    - FPlus   : flux vector from neighbour cell (N-side)
    //    - FMinus  : flux vector from owner cell (P-side)
    //    - n       : face unit normal (points from P to N)
    //    - UR      : conserved value on N-side at this Gauss point
    //    - UL      : conserved value on P-side at this Gauss point
    //  Output:
    //    - flux    : numerical flux vector (Cartesian)
    virtual void computeFlux
    (
        const label gaussID,
        const vector& FPlus,
        const vector& FMinus,
        const scalar UR,
        const scalar UL,
        const vector& n,
        vector& flux
    ) const override;

    //- Reading method
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //