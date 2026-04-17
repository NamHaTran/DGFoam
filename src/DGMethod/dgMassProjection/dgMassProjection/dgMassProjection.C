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

#include "dgMassProjection.H"
#include "PstreamReduceOps.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Helpers  * * * * * * * * * * * * * * * //

label dgMassProjection::findNameIndex
(
    const List<word>& names,
    const word& fieldName
)
{
    forAll(names, i)
    {
        if (names[i] == fieldName)
        {
            return i;
        }
    }

    return -1;
}


void dgMassProjection::validateCellID(const label cellID) const
{
    if (cellID < 0 || cellID >= mesh_.nCells())
    {
        FatalErrorInFunction
            << "Cell index " << cellID
            << " is out of range [0," << mesh_.nCells() - 1 << "]."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

dgMassProjection::dgMassProjection(const dgGeomMesh& mesh)
:
    mesh_(mesh),
    scalarNames_(),
    scalarFields_(),
    scalarResiduals_(),
    vectorNames_(),
    vectorFields_(),
    vectorResiduals_(),
    tensorNames_(),
    tensorFields_(),
    tensorResiduals_()
{}


// * * * * * * * * * * * * * Field Access  * * * * * * * * * * * * * * * * //

dgField<scalar>& dgMassProjection::scalarField
(
    const word& fieldName
)
{
    const label fieldI = findNameIndex(scalarNames_, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgMassProjection."
            << abort(FatalError);
    }

    return *scalarFields_[fieldI];
}


const dgField<scalar>& dgMassProjection::scalarField
(
    const word& fieldName
) const
{
    const label fieldI = findNameIndex(scalarNames_, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Scalar field " << fieldName
            << " is not registered in dgMassProjection."
            << abort(FatalError);
    }

    return *scalarFields_[fieldI];
}


dgField<vector>& dgMassProjection::vectorField
(
    const word& fieldName
)
{
    const label fieldI = findNameIndex(vectorNames_, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgMassProjection."
            << abort(FatalError);
    }

    return *vectorFields_[fieldI];
}


const dgField<vector>& dgMassProjection::vectorField
(
    const word& fieldName
) const
{
    const label fieldI = findNameIndex(vectorNames_, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Vector field " << fieldName
            << " is not registered in dgMassProjection."
            << abort(FatalError);
    }

    return *vectorFields_[fieldI];
}


dgField<tensor>& dgMassProjection::tensorField
(
    const word& fieldName
)
{
    const label fieldI = findNameIndex(tensorNames_, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Tensor field " << fieldName
            << " is not registered in dgMassProjection."
            << abort(FatalError);
    }

    return *tensorFields_[fieldI];
}


const dgField<tensor>& dgMassProjection::tensorField
(
    const word& fieldName
) const
{
    const label fieldI = findNameIndex(tensorNames_, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Tensor field " << fieldName
            << " is not registered in dgMassProjection."
            << abort(FatalError);
    }

    return *tensorFields_[fieldI];
}


List<List<scalar>>& dgMassProjection::scalarResidual
(
    const word& fieldName
)
{
    return residualBuffer<scalar>(fieldName);
}


const List<List<scalar>>& dgMassProjection::scalarResidual
(
    const word& fieldName
) const
{
    return residualBuffer<scalar>(fieldName);
}


List<List<vector>>& dgMassProjection::vectorResidual
(
    const word& fieldName
)
{
    return residualBuffer<vector>(fieldName);
}


const List<List<vector>>& dgMassProjection::vectorResidual
(
    const word& fieldName
) const
{
    return residualBuffer<vector>(fieldName);
}


List<List<tensor>>& dgMassProjection::tensorResidual
(
    const word& fieldName
)
{
    return residualBuffer<tensor>(fieldName);
}


const List<List<tensor>>& dgMassProjection::tensorResidual
(
    const word& fieldName
) const
{
    return residualBuffer<tensor>(fieldName);
}


// * * * * * * * * * * * * * Driver API  * * * * * * * * * * * * * * * * * //

void dgMassProjection::clearResiduals()
{
    forAll(scalarResiduals_, fieldI)
    {
        clearResidual(scalarResiduals_[fieldI]);
    }

    forAll(vectorResiduals_, fieldI)
    {
        clearResidual(vectorResiduals_[fieldI]);
    }

    forAll(tensorResiduals_, fieldI)
    {
        clearResidual(tensorResiduals_[fieldI]);
    }
}


void dgMassProjection::clearResiduals(const label cellID)
{
    validateCellID(cellID);

    forAll(scalarResiduals_, fieldI)
    {
        clearResidualCell(scalarResiduals_[fieldI], cellID);
    }

    forAll(vectorResiduals_, fieldI)
    {
        clearResidualCell(vectorResiduals_[fieldI], cellID);
    }

    forAll(tensorResiduals_, fieldI)
    {
        clearResidualCell(tensorResiduals_[fieldI], cellID);
    }
}


void dgMassProjection::clearResiduals(const word& fieldName)
{
    const label scalarI = findNameIndex(scalarNames_, fieldName);
    if (scalarI >= 0)
    {
        clearResidual(scalarResiduals_[scalarI]);
        return;
    }

    const label vectorI = findNameIndex(vectorNames_, fieldName);
    if (vectorI >= 0)
    {
        clearResidual(vectorResiduals_[vectorI]);
        return;
    }

    const label tensorI = findNameIndex(tensorNames_, fieldName);
    if (tensorI >= 0)
    {
        clearResidual(tensorResiduals_[tensorI]);
        return;
    }

    FatalErrorInFunction
        << "Field " << fieldName
        << " is not registered in dgMassProjection."
        << abort(FatalError);
}


void dgMassProjection::clearResiduals
(
    const word& fieldName,
    const label cellID
)
{
    validateCellID(cellID);

    const label scalarI = findNameIndex(scalarNames_, fieldName);
    if (scalarI >= 0)
    {
        clearResidualCell(scalarResiduals_[scalarI], cellID);
        return;
    }

    const label vectorI = findNameIndex(vectorNames_, fieldName);
    if (vectorI >= 0)
    {
        clearResidualCell(vectorResiduals_[vectorI], cellID);
        return;
    }

    const label tensorI = findNameIndex(tensorNames_, fieldName);
    if (tensorI >= 0)
    {
        clearResidualCell(tensorResiduals_[tensorI], cellID);
        return;
    }

    FatalErrorInFunction
        << "Field " << fieldName
        << " is not registered in dgMassProjection."
        << abort(FatalError);
}


void dgMassProjection::solve()
{
    solveType<scalar>();
    solveType<vector>();
    solveType<tensor>();
}


void dgMassProjection::solve(const label cellID)
{
    validateCellID(cellID);

    solveType<scalar>(cellID);
    solveType<vector>(cellID);
    solveType<tensor>(cellID);
}


void dgMassProjection::solve(const word& fieldName)
{
    const label scalarI = findNameIndex(scalarNames_, fieldName);
    if (scalarI >= 0)
    {
        solveField(*scalarFields_[scalarI], scalarResiduals_[scalarI]);
        return;
    }

    const label vectorI = findNameIndex(vectorNames_, fieldName);
    if (vectorI >= 0)
    {
        solveField(*vectorFields_[vectorI], vectorResiduals_[vectorI]);
        return;
    }

    const label tensorI = findNameIndex(tensorNames_, fieldName);
    if (tensorI >= 0)
    {
        solveField(*tensorFields_[tensorI], tensorResiduals_[tensorI]);
        return;
    }

    FatalErrorInFunction
        << "Field " << fieldName
        << " is not registered in dgMassProjection."
        << abort(FatalError);
}


void dgMassProjection::solve
(
    const word& fieldName,
    const label cellID
)
{
    validateCellID(cellID);

    const label scalarI = findNameIndex(scalarNames_, fieldName);
    if (scalarI >= 0)
    {
        solveFieldCell(*scalarFields_[scalarI], scalarResiduals_[scalarI], cellID);
        return;
    }

    const label vectorI = findNameIndex(vectorNames_, fieldName);
    if (vectorI >= 0)
    {
        solveFieldCell(*vectorFields_[vectorI], vectorResiduals_[vectorI], cellID);
        return;
    }

    const label tensorI = findNameIndex(tensorNames_, fieldName);
    if (tensorI >= 0)
    {
        solveFieldCell(*tensorFields_[tensorI], tensorResiduals_[tensorI], cellID);
        return;
    }

    FatalErrorInFunction
        << "Field " << fieldName
        << " is not registered in dgMassProjection."
        << abort(FatalError);
}


// * * * * * * * * * * * * * Diagnostics  * * * * * * * * * * * * * * * * * //

void dgMassProjection::writeInfo(Ostream& os) const
{
    os  << "dgMassProjection:" << nl
        << "  scalar fields : " << scalarNames_.size() << nl
        << "  vector fields : " << vectorNames_.size() << nl
        << "  tensor fields : " << tensorNames_.size() << nl;

    if (scalarNames_.size())
    {
        os << "    scalar : " << scalarNames_ << nl;
    }

    if (vectorNames_.size())
    {
        os << "    vector : " << vectorNames_ << nl;
    }

    if (tensorNames_.size())
    {
        os << "    tensor : " << tensorNames_ << nl;
    }

    os << nl;
}


void dgMassProjection::writeResiduals(Ostream& os) const
{
    os << "dgMassProjection residuals:" << nl;

    writeResidualType<scalar>(os, "scalar");
    writeResidualType<vector>(os, "vector");
    writeResidualType<tensor>(os, "tensor");

    os << nl;
}

} // End namespace Foam

// ************************************************************************* //
