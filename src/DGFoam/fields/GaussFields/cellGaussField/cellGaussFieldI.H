/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

#ifndef Foam_cellGaussFieldI_H
#define Foam_cellGaussFieldI_H

#include "cellGaussField.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * Constructors * * * * * * * * //

// Default
template<class Type>
cellGaussField<Type>::cellGaussField()
:
    cellID_(-1),
    dgMesh_(nullptr),
    cell_(nullptr),
    nGauss_(0),
    nDof_(0),
    dof_(nullptr),
    values_()
{}


// From DOF pointer
template<class Type>
cellGaussField<Type>::cellGaussField
(
    const dgGeomMesh* dgMesh,
    const cellDof<Type>* dof
)
:
    cellID_(dof ? dof->cellID() : -1),
    dgMesh_(dgMesh),
    cell_(dgMesh ? dgMesh->cells()[cellID_] : nullptr),
    nGauss_(cell_ ? cell_->gaussPoints().size() : 0),
    nDof_(dof ? dof->nDof() : 0),
    dof_(dof),
    values_(nGauss_)
{
    if (!dgMesh_)
    {
        FatalErrorInFunction << "dgMesh pointer is null" << abort(FatalError);
    }
}


// Copy constructor
template<class Type>
cellGaussField<Type>::cellGaussField
(
    const cellGaussField<Type>& other
)
:
    cellID_(other.cellID_),
    dgMesh_(other.dgMesh_),
    cell_(other.cell_),
    nGauss_(other.nGauss_),
    nDof_(other.nDof_),
    dof_(other.dof_),
    values_(other.values_)
{}


// Uniform value constructor
template<class Type>
cellGaussField<Type>::cellGaussField
(
    const label cellID,
    const dgGeomMesh* dgMesh,
    const Type& initVal
)
:
    cellID_(cellID),
    dgMesh_(dgMesh),
    cell_(dgMesh->cells()[cellID]),
    nGauss_(cell_->gaussPoints().size()),
    nDof_(0),
    dof_(nullptr),
    values_(nGauss_, initVal)
{
    if (!dgMesh_)
    {
        FatalErrorInFunction << "dgMesh pointer is null" << abort(FatalError);
    }
}


// Zero-init constructor
template<class Type>
cellGaussField<Type>::cellGaussField
(
    const label cellID,
    const dgGeomMesh* dgMesh
)
:
    cellID_(cellID),
    dgMesh_(dgMesh),
    cell_(dgMesh->cells()[cellID]),
    nGauss_(cell_->gaussPoints().size()),
    nDof_(0),
    dof_(nullptr),
    values_(nGauss_)
{
    if (!dgMesh_)
    {
        FatalErrorInFunction << "dgMesh pointer is null" << abort(FatalError);
    }
}


// * * * * * * * * Member Functions * * * * * * * * //

// Interpolate DOF → Gauss values
template<class Type>
void cellGaussField<Type>::interpolateFromDof()
{
    if (!dof_)
    {
        FatalErrorInFunction << "No DOF data to interpolate" << abort(FatalError);
    }

    const List<List<scalar>>& B = cell_->basis();

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        Type val = pTraits<Type>::zero;

        for (label k = 0; k < nDof_; ++k)
        {
            val += dof_->dof()[k] * B[gp][k];
        }

        values_[gp] = val;
    }
}


// Copy assignment
template<class Type>
inline cellGaussField<Type>& cellGaussField<Type>::operator=
(
    const cellGaussField<Type>& other
)
{
    if (this == &other) return *this;

    cellID_ = other.cellID_;
    dgMesh_ = other.dgMesh_;
    cell_   = other.cell_;
    nGauss_ = other.nGauss_;
    nDof_   = other.nDof_;
    dof_    = other.dof_;
    values_ = other.values_;

    return *this;
}


// Assign uniform value
template<class Type>
inline cellGaussField<Type>& cellGaussField<Type>::operator=
(
    const Type& value
)
{
    for (label gp = 0; gp < nGauss_; ++gp)
    {
        values_[gp] = value;
    }
    return *this;
}

} // namespace Foam

#endif

// ************************************************************************* //

