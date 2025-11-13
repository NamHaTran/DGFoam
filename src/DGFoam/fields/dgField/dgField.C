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

\*---------------------------------------------------------------------------*/

#include "dgField.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::dgField<Type>::dgField
(
    const IOobject& io,
    const dgGeomMesh& mesh,
    bool hasDof
)
:
    regIOobject(io),
    mesh_(mesh),
    hasDof_(hasDof),
    dofPtr_(nullptr),
    bcManagerPtr_(nullptr)
{
    if (hasDof)
    {
        // Create DOF field from IOobject
        dofPtr_ = new dofField<Type>
        (
            io.name(),
            io.instance(),
            io.readOpt(),
            io.writeOpt(),
            mesh
        );

        gaussFields_.setSize(mesh_.nCells());

        // Interpolate Gauss values for each cell
        forAll(gaussFields_, cellI)
        {
            gaussFields_[cellI] = GaussField<Type>
            (
                dofPtr_,
                cellI,
                &mesh
            );

            // Compute values at Gauss points
            gaussFields_[cellI].interpolateFromDof();
        }

        Info<< "dgField<" << pTraits<Type>::typeName
            << ">: constructed conservative field " << io.name()
            << " with " << mesh.nCells() << " cells." << nl;
    }
    else
    {
        // Create empty GaussFields only (no DOF)
        forAll(gaussFields_, cellI)
        {
            gaussFields_[cellI] = GaussField<Type>(nullptr, cellI, &mesh);
        }

        Info<< "dgField<" << pTraits<Type>::typeName
            << ">: constructed primary field " << io.name()
            << " (no DOF data)." << nl;
    }
}

// Destructor
template<class Type>
Foam::dgField<Type>::~dgField()
{
    // Only delete DOF field if this dgField owns it
    if (hasDof_ && dofPtr_)
    {
        delete dofPtr_;
        dofPtr_ = nullptr;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::dgField<Type>::correctBoundary()
{
    // Placeholder function
    // Will be implemented when boundary conditions are integrated
}


template<class Type>
bool Foam::dgField<Type>::writeData(Ostream& os) const
{
    /*
    if (dofPtr_.valid())
    {
        // Delegate writing to DOF field
        return dofPtr_().writeData(os);
    }

    // For primary fields, no DOF data to write
    Info<< "dgField<" << pTraits<Type>::typeName
        << ">: no DOF data to write for " << name() << nl;
    */

    return true;
}

// * * * * * * * * * * * * * * Template instantiation * * * * * * * * * * * * //
template class Foam::dgField<scalar>;
template class Foam::dgField<vector>;
template class Foam::dgField<tensor>;
template class Foam::dgField<symmTensor>;
template class Foam::dgField<sphericalTensor>;

// Explicit instantiation of static member typeName
template<> const word dgField<scalar>::typeName("dgField<scalar>");
template<> const word dgField<vector>::typeName("dgField<vector>");
template<> const word dgField<tensor>::typeName("dgField<tensor>");
template<> const word dgField<symmTensor>::typeName("dgField<symmTensor>");
template<> const word dgField<sphericalTensor>::typeName("dgField<sphericalTensor>");

// ************************************************************************* //

} // End namespace Foam

