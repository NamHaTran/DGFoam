/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::rhoBasedConservative

Description
    Conservative thermodynamics model for DG compressible solvers.

    Reconstructs all primitive thermo quantities (T, p, e, h, Cp, Cv, gamma,
    a, mu, kappa, Pr) at each Gauss point from the conservative fields
    rho, rhoU, and E. Supports ideal-gas and thermal-perfect-gas EOS and
    compatible thermo/transport/energy laws.

    Intended for DG solvers using the conservative Navier–Stokes formulation.
    Not applicable to primitive-variable or incompressible formulations.

SourceFiles
    rhoBasedConservative.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_rhoBasedConservative_H
#define Foam_rhoBasedConservative_H

#include "dgThermoConservative.H"
#include "eqnOfState.H"
#include "thermoLaw.H"
#include "transportLaw.H"
#include "energy.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class rhoBasedConservative Declaration
\*---------------------------------------------------------------------------*/

class rhoBasedConservative
:
    public dgThermoConservative
{
public:

    //- Runtime type information
    TypeName("rhoBasedConservative");

    //- Runtime selection
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoBasedConservative,
        dictionary,
        (const word& name, const dictionary& dict, const dgGeomMesh& mesh),
        (name, dict, mesh)
    );


    // Constructors

        //- Disallow default constructor
        rhoBasedConservative() = delete;

        //- Construct from name, dictionary and DG mesh
        rhoBasedConservative
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        );

    //- Destructor
    virtual ~rhoBasedConservative() = default;


    // Model setup interface

        //- Construct EOS, thermo, transport and energy models
        virtual void initModels() override;

        //- Validate model combination (EOS/thermo/transport/energy)
        virtual void validateModels() override;


    // Core API

        //- Update thermo state (cell-wise interface)
        //  - For now this forwards to correct(), more optimized
        //    per-cell version can be added later if needed
        virtual void update(const label& cellID) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
