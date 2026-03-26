/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgLimiterManager.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "error.H"

namespace Foam
{

IOobject dgLimiterManager::createIOobject(const dgGeomMesh& mesh) const
{
    IOobject io
    (
        "dgOptions",
        mesh.getFvMesh().time().constant(),
        mesh.getFvMesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        io.readOpt(IOobject::NO_READ);
    }

    return io;
}


void dgLimiterManager::constructLimiters()
{
    limiters_.clear();
    limiterNames_.clear();

    if (readOpt() == IOobject::NO_READ)
    {
        return;
    }

    if (!found("limiters"))
    {
        return;
    }

    const dictionary& limitersDict = subDict("limiters");
    DynamicList<word> limiterNames(limitersDict.size());

    forAllConstIter(dictionary, limitersDict, iter)
    {
        const word key = iter().keyword();

        if (key == "FoamFile")
        {
            continue;
        }

        if (!limitersDict.isDict(key))
        {
            FatalIOErrorInFunction(limitersDict)
                << "Limiter entry '" << key
                << "' must be a sub-dictionary."
                << exit(FatalIOError);
        }

        limiterNames.append(key);
    }

    limiterNames_.transfer(limiterNames);
    limiters_.setSize(limiterNames_.size());

    forAll(limiterNames_, limiterI)
    {
        const word& limiterName = limiterNames_[limiterI];
        autoPtr<dgLimiter> limiter =
            dgLimiter::New(limitersDict.subDict(limiterName), mesh_);

        limiters_.set(limiterI, limiter.ptr());
    }
}


dgLimiterManager::dgLimiterManager
(
    const dgGeomMesh& mesh
)
:
    IOdictionary(createIOobject(mesh)),
    mesh_(mesh),
    limiters_(),
    limiterNames_()
{
    constructLimiters();
}


bool dgLimiterManager::has(const word& name) const
{
    forAll(limiterNames_, limiterI)
    {
        if (limiterNames_[limiterI] == name)
        {
            return true;
        }
    }

    return false;
}


const dgLimiter& dgLimiterManager::limiter(const word& name) const
{
    forAll(limiterNames_, limiterI)
    {
        if (limiterNames_[limiterI] == name)
        {
            return limiters_[limiterI];
        }
    }

    FatalErrorInFunction
        << "Unknown limiter name '" << name << "'." << nl
        << "Configured limiter names: " << limiterNames_
        << exit(FatalError);

    return limiters_[0];
}


dgLimiter& dgLimiterManager::limiter(const word& name)
{
    forAll(limiterNames_, limiterI)
    {
        if (limiterNames_[limiterI] == name)
        {
            return limiters_[limiterI];
        }
    }

    FatalErrorInFunction
        << "Unknown limiter name '" << name << "'." << nl
        << "Configured limiter names: " << limiterNames_
        << exit(FatalError);

    return limiters_[0];
}


void dgLimiterManager::correct()
{
    forAll(limiters_, limiterI)
    {
        limiters_[limiterI].correct();
    }
}


bool dgLimiterManager::read()
{
    if (IOdictionary::regIOobject::read())
    {
        constructLimiters();
        return true;
    }

    return false;
}


void dgLimiterManager::listLimiters(Ostream& os) const
{
    os  << "dgLimiterManager: configured limiters = "
        << limiterNames_.size() << nl;

    forAll(limiterNames_, limiterI)
    {
        word limiterType("unknown");
        const dictionary& limiterDict = limiters_[limiterI].dict();

        if (limiterDict.found("type"))
        {
            limiterType = limiterDict.get<word>("type");
        }

        os  << "  - " << limiterNames_[limiterI]
            << " : " << limiterType << nl;
    }

    os << nl;
}

} // End namespace Foam

// ************************************************************************* //
