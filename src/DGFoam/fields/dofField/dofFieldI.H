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
    Foam::dofField<Type>

Description
    Store degrees of freedom (Dof) of given Type for a DG cell.

SourceFiles
    dofFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dofFieldI_H
#define Foam_dofFieldI_H

#include "dofField.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * * * Constructors * * * * * * * * * * * * //

template<class Type>
dofField<Type>::dofField
(
    const word& name,
    const fileName& instance,
    IOobject::readOption rOpt,
    IOobject::writeOption wOpt,
    const dgGeomMesh& mesh
)
:
    dgMesh_(mesh),
    cellDofs_(),
    nCells_(mesh.nCells()),
    foamFields_(),
    dofPerCell_(nCells_, 0),
    maxDoF_(0)
{
    // Determine number of DoFs per cell and maxDoF_
    for (label c = 0; c < nCells_; ++c)
    {
        const label nDof = mesh.cells()[c]->nDof();
        dofPerCell_[c] = nDof;
        maxDoF_ = max(maxDoF_, nDof);
    }

    // Read volFields (one per DoF index)
    foamFields_.setSize(maxDoF_);

    for (label d = 0; d < maxDoF_; ++d)
    {
        const word fieldName = name + "_" + Foam::name(d);

        IOobject probeIO
        (
            fieldName,
            instance,
            mesh.getFvMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        // Determine if volField file exists
        if (probeIO.typeHeaderOk<GeometricField<Type, fvPatchField, volMesh>>(true))
        {
            Info<< "Reading DoF volField: " << fieldName << nl;

            foamFields_[d].reset
            (
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fieldName,
                        instance,
                        mesh.getFvMesh(),
                        IOobject::MUST_READ,
                        wOpt
                    ),
                    mesh.getFvMesh()
                )
            );
        }
        else
        {
            Info<< "DoF volField not found: " << fieldName
                << " → creating zero field." << nl;

            foamFields_[d].reset
            (
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fieldName,
                        instance,
                        mesh.getFvMesh(),
                        IOobject::NO_READ,
                        wOpt
                    ),
                    mesh.getFvMesh(),
                    dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
                )
            );
        }
    }

    // Construct cellDof list from volFields
    cellDofs_.setSize(nCells_);

    for (label c = 0; c < nCells_; ++c)
    {
        const label nDof = dofPerCell_[c];
        List<Type> local(nDof, pTraits<Type>::zero);

        for (label d = 0; d < nDof; ++d)
        {
            local[d] = (*foamFields_[d])[c];
        }

        cellDofs_[c] = cellDof<Type>(c, nDof, local);
    }
}


// * * * * * * * * * * Member Functions  * * * * * * * * * //

template<class Type>
void dofField<Type>::updateFoamFields()
{
    for (label c = 0; c < nCells_; ++c)
    {
        const label nDof = dofPerCell_[c];

        for (label d = 0; d < maxDoF_; ++d)
        {
            if (d < nDof)
            {
                (*foamFields_[d])[c] = cellDofs_[c][d];
            }
            else
            {
                (*foamFields_[d])[c] = pTraits<Type>::zero;
            }
        }
    }
}


template<class Type>
void dofField<Type>::write() const
{
    for (const auto& fld : foamFields_)
    {
        fld->write();
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
