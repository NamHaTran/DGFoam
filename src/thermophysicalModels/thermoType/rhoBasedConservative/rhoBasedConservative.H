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
    DG pipeline: density-based conservative
    - Inputs: conserved (rhoC, rhoU, rhoE)
    - Recover primitives via EOS + thermoLaw
    - Compute transport via transportLaw
    - Restrict inner models to:
        * eqnOfState: idealGas
        * thermoLaw     : hConst (constant Cp)
        * transportLaw  : Sutherland, powerVHS  (extendable)
        * energy        : sensibleInternalEnergy

SourceFiles
    rhoBasedConservative.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_rhoBasedConservativeThermo_H
#define Foam_rhoBasedConservativeThermo_H

#include "dgThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "GaussField.H"

namespace Foam
{

class rhoBasedConservative : public dgThermo
{
private:

    // Attributes

        //- Gas constant (cached from eqnState)
        scalar Rgas_;

public:

    //- Runtime type information
    TypeName("rhoBasedConservative");

    //- Runtime selection
    declareRunTimeSelectionTable
    (
        autoPtr, rhoBasedConservative, dictionary,
        (const word& name, const dictionary& dict, dgGeomMesh& mesh),
        (name, dict, mesh)
    );

    //- Construct
    rhoBasedConservative(const word& name, const dictionary& dict, dgGeomMesh& mesh);

    //- Destructor
    ~rhoBasedConservative() override = default;

    //- Build inner models
    void initModels() override;

    //- Validate allowed inner model combination
    void validateModelCombo() override;

    //- Compute properties from conserved inputs
    void update
    (
        const dgThermoInputs& in,
        dgThermoOutputs& out
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
