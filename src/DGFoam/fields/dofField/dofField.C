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

#include "dofField.H"
#include "volFields.H"

namespace Foam
{

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
    // Determine dof per cell and maxDoF_
    for (label c = 0; c < nCells_; ++c)
    {
        label nDof = mesh.cells()[c]->nDof();
        dofPerCell_[c] = nDof;
        maxDoF_ = max(maxDoF_, nDof);
    }

    // Read foam fields from disk
    foamFields_.setSize(maxDoF_);
    for (label d = 0; d < maxDoF_; ++d)
    {
        word fieldName = name + "_" + Foam::name(d);

        IOobject io
        (
            fieldName,
            instance,
            mesh.getFvMesh(),
            IOobject::READ_IF_PRESENT,  // Temporary: just to check availability
            IOobject::NO_WRITE
        );

        if (io.typeHeaderOk<GeometricField<Type, fvPatchField, volMesh>>(true))  // File exists and is readable
        {
            Info<< "Reading field from file: " << fieldName << nl;

            foamFields_[d].reset
            (
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fieldName,
                        instance,
                        mesh.getFvMesh(),
                        IOobject::MUST_READ,   // Force read
                        wOpt
                    ),
                    mesh.getFvMesh()
                )
            );
        }
        else
        {
            Info<< "File not found for field " << fieldName << ". Creating a new zero field instead." << nl;

            foamFields_[d].reset
            (
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fieldName,
                        instance,
                        mesh.getFvMesh(),
                        IOobject::NO_READ,     // No read
                        wOpt                   // Still allow writing later
                    ),
                    mesh.getFvMesh(),
                    dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
                )
            );
        }
    }


    // Construct cellDofs from foamFields
    cellDofs_.setSize(nCells_);
    for (label c = 0; c < nCells_; ++c)
    {
        label nDof = dofPerCell_[c];

        List<Type> inputDof(nDof, pTraits<Type>::zero);

        for (label d = 0; d < nDof; ++d)
        {
            inputDof[d] = (*foamFields_[d])[c];
        }

        cellDofs_[c] = cellDof<Type>(c, nDof, inputDof);
    }
}


template<class Type>
void dofField<Type>::updateFoamFields()
{
    for (label c = 0; c < nCells_; ++c)
    {
        label nDof = dofPerCell_[c];

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

// * * * * * * * * * * * * * * Template instantiation * * * * * * * * * * * * //
template class Foam::dofField<scalar>;
template class Foam::dofField<vector>;
template class Foam::dofField<tensor>;
template class Foam::dofField<symmTensor>;
template class Foam::dofField<sphericalTensor>;

} // End namespace Foam

// ************************************************************************* //
