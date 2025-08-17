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

#include "dgThermo.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dgThermo, 0);
defineRunTimeSelectionTable(dgThermo, dictionary);

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

Foam::dgThermo::dgThermo
(
    const word& name,
    const dictionary& dict,
    dgGeomMesh& mesh
)
:
    name_(name),
    dict_(dict),
    mesh_(mesh),
    eqnState_(nullptr),
    thermo_(nullptr),
    transport_(nullptr),
    energy_(nullptr),
    R_(Zero), Cp_(Zero), h_(Zero), e_(Zero),
    mu_(Zero), kappa_(Zero), Pr_(Zero), a_(Zero)
{
    // Base: no model construction here (defer to derived::initModels()).
}

// * * * * * * * * * * * * * Factory Method  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dgThermo> Foam::dgThermo::New
(
    const word& name,
    const dictionary& dict,
    dgGeomMesh& mesh
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgThermo type: " << name << nl
            << "Valid dgThermo types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dict, mesh);
}

// ************************************************************************* //

} // End namespace Foam
