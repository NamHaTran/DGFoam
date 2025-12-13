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

Description
    Central flux for DG:
        h = 0.5 * (F_plus + F_minus)

    No wave-speed estimate, no dissipation.
    Suitable for diffusive fluxes or LDG auxiliary equations.

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

class dgCentralFluxSolver : public dgFluxSolver
{
public:

    //- Runtime type name
    TypeName("central");


    // Constructors
    dgCentralFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    )
    :
        dgFluxSolver(name, dict, mesh)
    {}


    virtual ~dgCentralFluxSolver() = default;


    //- For scalar U and vector F
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;


    //- For vector U and tensor F
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;


    //- No state to reset
    virtual void reset() override {}

    //- No settings to read
    virtual void read(const dictionary& dict) override {}
};

} // End namespace Foam

#endif

// ************************************************************************* //

