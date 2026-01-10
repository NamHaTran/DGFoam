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
    Foam::dgHLLCFluxSolver

Description
    Local Lax–Friedrichs / Rusanov numerical flux.

    Numerical flux:
        h(U+, U-, n) =
            0.5*(F(U+) + F(U-))·n
          - 0.5*C*(U+ - U-)

    with dissipation coefficient:
        C = max( |U_n^-| + a^- , |U_n^+| + a^+ )

    Optional:
        If 'scaleByMach = true' in dictionary:
            C = |U_n| + a/M     where  M = |U_n|/a

    Notes:
    - For scalar transport: Cartesian flux = fn * n
    - For momentum transport: output flux is a vector, not multiplied by n

SourceFiles
    dgHLLCFluxSolver.C
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

/*---------------------------------------------------------------------------*\
                Class dgLaxFriedrichsFluxSolver Declaration
\*---------------------------------------------------------------------------*/

namespace Foam
{

class dgLaxFriedrichsFluxSolver : public dgFluxSolver
{
private:

    //- Velocity field
    const dgField<vector>& U_;

    //- Speed of sound field
    const dgField<scalar>& a_;

    //- Enable Mach-based dissipation model
    bool scaleByMach_;


    //- Compute dissipation coefficient C
    scalar calcDissipationCoeff
    (
        const vector& ULv,
        const vector& URv,
        const scalar aL,
        const scalar aR,
        const vector& n
    ) const;

public:

    //- Runtime type name
    TypeName("LaxFriedrichs");

    //- Construct from name, dictionary and mesh
    dgLaxFriedrichsFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    virtual ~dgLaxFriedrichsFluxSolver() = default;


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

    //- No caching required
    virtual void reset() override {}

    //- Read dictionary
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
