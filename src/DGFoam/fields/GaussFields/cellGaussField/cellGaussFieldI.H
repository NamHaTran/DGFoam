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

// Move constructor
template<class Type>
cellGaussField<Type>::cellGaussField
(
    cellGaussField<Type>&& other
) noexcept
:
    cellID_(other.cellID_),
    dgMesh_(other.dgMesh_),
    cell_(other.cell_),
    nGauss_(other.nGauss_),
    nDof_(other.nDof_),
    dof_(other.dof_),
    values_(std::move(other.values_))
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

/*
template<class Type>
cellGaussField<Type>::cellGaussField(const tmp<cellGaussField<Type>>& t)
{
    if (t.isTmp())
    {
        this->operator=(std::move(const_cast<cellGaussField<Type>&>(t())));
        t.clear();
    }
    else
    {
        *this = t();
    }
}

template<class Type>
Foam::cellGaussField<Type>::cellGaussField
(
    const tmp<cellGaussField<Type>>& t
)
:
    cellGaussField<Type>(t())
{
    // Release ownership from tmp to avoid deletion of the internal object
    t.clear();
}
*/

template<class Type>
Foam::cellGaussField<Type>::cellGaussField
(
    const tmp<cellGaussField<Type>>& t
)
{
    cellGaussField<Type>& tf =
        const_cast<tmp<cellGaussField<Type>>&>(t).ref();

    cellID_ = tf.cellID_;
    dgMesh_ = tf.dgMesh_;
    cell_   = tf.cell_;
    nGauss_ = tf.nGauss_;
    nDof_   = tf.nDof_;
    dof_    = tf.dof_;

    // Transfer storage instead of copy
    values_.transfer(tf.values_);

    t.clear();
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

// Move assignment
template<class Type>
inline cellGaussField<Type>& cellGaussField<Type>::operator=
(
    cellGaussField<Type>&& other
) noexcept
{
    if (this != &other)
    {
        cellID_ = other.cellID_;
        dgMesh_ = other.dgMesh_;
        cell_   = other.cell_;
        nGauss_ = other.nGauss_;
        nDof_   = other.nDof_;
        dof_    = other.dof_;

        values_ = std::move(other.values_);
    }
    return *this;
}

// Assign tmp
template<class Type>
cellGaussField<Type>& Foam::cellGaussField<Type>::operator=
(
    const tmp<cellGaussField<Type>>& t
)
{
    /*
    *this = t();
    t.clear();
    return *this;
    */
   
    cellGaussField<Type>& tf =
        const_cast<tmp<cellGaussField<Type>>&>(t).ref();

    cellID_ = tf.cellID_;
    dgMesh_ = tf.dgMesh_;
    cell_   = tf.cell_;
    nGauss_ = tf.nGauss_;
    nDof_   = tf.nDof_;
    dof_    = tf.dof_;

    values_.transfer(tf.values_);

    t.clear();

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

