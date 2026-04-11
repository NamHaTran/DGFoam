/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Inline implementations for Foam::dgField.

SourceFiles
    dgFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFieldI_H
#define Foam_dgFieldI_H

#include "dgField.H"
#include "IOstreams.H"
#include "Jacobian.H"
#include "dgProcessorBoundaryManager.H"

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
    dofPtr_(nullptr),
    procSyncPtr_(nullptr)
{
    const label nCells = mesh_.nCells();
    gaussFields_.setSize(nCells);

    if (hasDof)
    {
        // Create DOF container
        dofPtr_.reset
        (
            new dofField<Type>
            (
                io.name(),
                io.instance(),
                io.readOpt(),
                io.writeOpt(),
                mesh
            )
        );

        // Interpolate DOF → Gauss values
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            gaussFields_[cellI] = GaussField<Type>(dofPtr_, cellI, &mesh_);
            gaussFields_[cellI].interpolateFromDof();
        }

        //Info<< "dgField<" << pTraits<Type>::typeName
        //    << ">: constructed conservative field " << io.name()
        //    << " with " << nCells << " cells." << nl;
    }
    else
    {
        // Primary field: allocate Gauss structure only
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            gaussFields_[cellI] = GaussField<Type>(cellI, &mesh_);
        }

        //Info<< "dgField<" << pTraits<Type>::typeName
        //    << ">: constructed primary field " << io.name()
        //    << " (no DOF)." << nl;
    }
}

// * * * * * * * * * Member Functions  * * * * * * * * //

template<class Type>
void dgField<Type>::synch()
{
    if (!procSyncPtr_.valid())
    {
        procSyncPtr_.reset(new dgProcessorBoundaryManager<Type>(mesh_));
    }

    procSyncPtr_().update(*this);
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
dgField<Type>::cellAverageVolField() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tAvgField
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                this->name(),
                this->instance(),
                mesh_.getFvMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_.getFvMesh(),
            dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
        )
    );
    auto& avgField = tAvgField.ref();

    forAll(gaussFields_, cellI)
    {
        const dgGeomCell& cell = *mesh_.cells()[cellI];
        const cellGaussField<Type>& cellField = gaussFields_[cellI].cellField();
        const List<scalar>& weights = cell.weights();
        const List<scalar>& J3D = cell.J3D();

        Type integral = pTraits<Type>::zero;
        scalar totalWeight = 0.0;

        forAll(weights, gpI)
        {
            // J3D is det(dX/deta), including collapsed-coordinate factors.
            const scalar quadWeight = weights[gpI]*J3D[gpI];
            integral += cellField[gpI]*quadWeight;
            totalWeight += quadWeight;
        }

        if (mag(totalWeight) <= SMALL)
        {
            FatalErrorInFunction
                << "Non-positive quadrature weight while writing averaged field "
                << this->name() << " in cell " << cellI << nl
                << abort(FatalError);
        }

        avgField[cellI] = integral/totalWeight;
    }

    return tAvgField;
}

template<class Type>
bool dgField<Type>::writeData(Ostream& os) const
{
    return cellAverageVolField().cref().writeData(os);
}

template<class Type>
bool dgField<Type>::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if (dofPtr_)
    {
        dofPtr_->write();
    }

    tmp<GeometricField<Type, fvPatchField, volMesh>> tAvgField =
        cellAverageVolField();

    typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
        boundaryField = tAvgField.ref().boundaryFieldRef();

    forAll(boundaryField, patchI)
    {
        boundaryField[patchI] = boundaryField[patchI].patchInternalField();
    }

    return tAvgField.ref().writeObject(streamOpt, writeOnProc);
}

template<class Type>
dgField<Type>& dgField<Type>::operator=
(
    const Type& value
)
{
    const label nCells = mesh_.nCells();
    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        gaussFields_[cellI] = value;
    }
    return *this;
}

} // End namespace Foam

#endif
// ************************************************************************* //
