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

\*---------------------------------------------------------------------------*/

#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dgThermoConservative, 0);
defineRunTimeSelectionTable(dgThermoConservative, dictionary);

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

Foam::dgThermoConservative::dgThermoConservative
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    name_(name),
    dict_(dict),
    mesh_(mesh),

    // Null initialization. Derived class will allocate real models.
    eqnState_(nullptr),
    thermo_(nullptr),
    transport_(nullptr),
    energy_(nullptr),
    R_(Zero),
    Cp_
    (
        IOobject
        (
            "Cp", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    Cv_
    (
        IOobject
        (
            "Cv", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    h_
    (
        IOobject
        (
            "h", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    e_
    (
        IOobject
        (
            "e", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    mu_
    (
        IOobject
        (
            "mu", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    kappa_
    (
        IOobject
        (
            "kappa", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    Pr_
    (
        IOobject
        (
            "Pr", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    a_
    (
        IOobject
        (
            "a", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    gamma_
    (
        IOobject
        (
            "gamma", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    rho_
    (
        mesh_.getFvMesh().lookupObject<dgField<scalar>>("rho")
    ),
    rhoU_
    (
        mesh_.getFvMesh().lookupObject<dgField<vector>>("rhoU")
    ),
    E_
    (
        mesh_.getFvMesh().lookupObject<dgField<scalar>>("E")
    ),
    p_
    (
        mesh_.getFvMesh().lookupObjectRef<dgField<scalar>>("p")
    ),
    T_
    (
        mesh_.getFvMesh().lookupObjectRef<dgField<scalar>>("T")
    )
{
    // No model construction here.
    // Derived classes will call initModels() after base constructor completes.
}


// * * * * * * * * * * *  Factory Selector  * * * * * * * * * * * * //

Foam::autoPtr<Foam::dgThermoConservative> Foam::dgThermoConservative::New
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgThermoConservative type: " << name << nl
            << "Valid dgThermoConservative types: " << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dict, mesh);
}

} // End namespace Foam

// ************************************************************************* //

