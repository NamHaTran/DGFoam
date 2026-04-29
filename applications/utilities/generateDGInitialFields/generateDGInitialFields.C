/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Generate pOrder-0 DG conservative initial fields rho/rhoU/E from the
    uniform primitive internal state p/T/U in constant/boundaryConditions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "dgField.H"
#include "dgGeomMesh.H"
#include "energy.H"
#include "eqnOfState.H"
#include "thermoLaw.H"

namespace Foam
{
namespace
{

template<class Type>
void setCellDof0
(
    dgField<Type>& field,
    const label cellI,
    const Type& value
)
{
    cellDof<Type>& dofs = field.dof()[cellI];

    if (dofs.nDof() == 0)
    {
        FatalErrorInFunction
            << "Field " << field.name() << " has no DoFs in cell "
            << cellI << exit(FatalError);
    }

    dofs[0] = value;

    for (label dofI = 1; dofI < dofs.nDof(); ++dofI)
    {
        dofs[dofI] = pTraits<Type>::zero;
    }

    field.dof().updateCellDof(cellI);
    field.gaussFields()[cellI].interpolateFromDof();
}

template<class Type>
void writeDof0Field(const dgField<Type>& field, const dgGeomMesh& dgMesh)
{
    const word fieldName = field.name() + "_dof0";

    GeometricField<Type, fvPatchField, volMesh> dof0Field
    (
        IOobject
        (
            fieldName,
            field.instance(),
            dgMesh.getFvMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dgMesh.getFvMesh(),
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
    );

    for (label cellI = 0; cellI < field.nCells(); ++cellI)
    {
        dof0Field[cellI] = field.dof()[cellI][0];
    }

    typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
        boundaryField = dof0Field.boundaryFieldRef();

    forAll(boundaryField, patchI)
    {
        boundaryField[patchI] = boundaryField[patchI].patchInternalField();
    }

    Info<< "Writing " << fieldName << " to " << field.instance() << nl;
    dof0Field.write();
}

scalar calcStoredThermoEnergy
(
    const scalar rho,
    const scalar T,
    const eqnOfState& eos,
    const energy& energyModel
)
{
    if (energyModel.heIsInternalEnergy() && eos.canCalcEFromRhoT())
    {
        return eos.calcEFromRhoT(rho, T);
    }

    return energyModel.calcHe(T);
}

} // End namespace
} // End namespace Foam

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate DG conservative dof0 fields rho/rhoU/E at time 0 "
        "from constant/boundaryConditions internal p/T/U."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    runTime.setTime(instant(0, "0"), 0);

    Info<< "Create DG geometric mesh\n" << endl;
    dgGeomMesh dgMesh(mesh);

    Info<< "Create conservative DG fields\n" << endl;

    dgField<scalar> rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dgMesh,
        true
    );

    dgField<vector> rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dgMesh,
        true
    );

    dgField<scalar> E
    (
        IOobject
        (
            "E",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dgMesh,
        true
    );

    Info<< "Read thermophysical properties\n" << endl;

    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary thermoType(thermoDict.subDict("dgThermo"));
    const dictionary& mixDict = thermoDict.subDict("mixture");

    const word eosType(thermoType.lookup("equationOfState"));
    const word thermoLawType(thermoType.lookup("thermo"));
    const word energyType(thermoType.lookup("energy"));

    autoPtr<eqnOfState> eos =
        eqnOfState::New(eosType, mixDict, dgMesh);

    autoPtr<thermoLaw> thermo =
        thermoLaw::New(thermoLawType, mixDict, dgMesh, eos());

    autoPtr<energy> energyModel =
        energy::New(energyType, mixDict, dgMesh, thermo());

    Info<< "Read DG compressible boundary initial state\n" << endl;

    IOdictionary boundaryConditionsDict
    (
        IOobject
        (
            "boundaryConditions",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& internalDict =
        boundaryConditionsDict.subDict("internalField");

    const scalar p0 = internalDict.get<scalar>("p");
    const scalar T0 = internalDict.get<scalar>("T");
    const vector U0 = internalDict.get<vector>("U");

    Info<< "Initial primitive state: p = " << p0
        << ", T = " << T0
        << ", U = " << U0 << nl << endl;

    const scalar rho0 = eos->calcRhoFromPT(p0, T0);
    const scalar he0 = calcStoredThermoEnergy(rho0, T0, eos(), energyModel());
    const vector rhoU0 = rho0*U0;
    const scalar E0 = rho0*(he0 + 0.5*magSqr(U0));

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        setCellDof0(rho, cellI, rho0);
        setCellDof0(rhoU, cellI, rhoU0);
        setCellDof0(E, cellI, E0);
    }

    writeDof0Field(rho, dgMesh);
    writeDof0Field(rhoU, dgMesh);
    writeDof0Field(E, dgMesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
