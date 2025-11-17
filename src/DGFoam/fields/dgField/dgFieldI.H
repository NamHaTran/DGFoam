/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::dgField

Description
    DG field class storing DOF and GaussField data for the whole mesh.

SourceFiles
    dgBasisFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFieldI_H
#define Foam_dgFieldI_H

#include "dgField.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * Constructors  * * * * * * * * * //

template<class Type>
dgField<Type>::dgField
(
    const IOobject& io,
    const dgGeomMesh& mesh,
    bool hasDof
)
:
    regIOobject(io),
    mesh_(mesh),
    hasDof_(hasDof),
    dofPtr_(nullptr)
{
    const label nCells = mesh_.nCells();
    gaussFields_.setSize(nCells);

    if (hasDof)
    {
        // Create DOF container
        dofPtr_ = new dofField<Type>
        (
            io.name(),
            io.instance(),
            io.readOpt(),
            io.writeOpt(),
            mesh
        );

        // Interpolate DOF → Gauss values
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            gaussFields_[cellI] = GaussField<Type>(dofPtr_, cellI, &mesh_);
            gaussFields_[cellI].interpolateFromDof();
        }

        Info<< "dgField<" << pTraits<Type>::typeName
            << ">: constructed conservative field " << io.name()
            << " with " << nCells << " cells." << nl;
    }
    else
    {
        // Primary field: allocate Gauss structure only
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            gaussFields_[cellI] = GaussField<Type>(cellI, &mesh_);
        }

        Info<< "dgField<" << pTraits<Type>::typeName
            << ">: constructed primary field " << io.name()
            << " (no DOF)." << nl;
    }
}


// * * * * * * * * * Destructor * * * * * * * * * //

template<class Type>
dgField<Type>::~dgField()
{
    if (hasDof_ && dofPtr_)
    {
        delete dofPtr_;
        dofPtr_ = nullptr;
    }
}


// * * * * * * * * * Member Functions  * * * * * * * * //

template<class Type>
bool dgField<Type>::writeData(Ostream& os) const
{
    // For DOF field, its own write() handles write
    return true;
}

} // End namespace Foam

#endif
// ************************************************************************* //

