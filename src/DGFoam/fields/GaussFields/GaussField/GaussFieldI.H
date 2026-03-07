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

#ifndef Foam_GaussFieldI_H
#define Foam_GaussFieldI_H

#include "GaussField.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * Constructors * * * * * * * * //

template<class Type>
GaussField<Type>::GaussField()
:
    dofField_(nullptr),
    cellID_(-1),
    mesh_(nullptr),
    cellField_(),
    faceField_()
{}


template<class Type>
GaussField<Type>::GaussField
(
    const dofField<Type>* dofField,
    label cellID,
    const dgGeomMesh* mesh
)
:
    dofField_(dofField),
    cellID_(cellID),
    mesh_(mesh),
    cellField_(mesh, &(*dofField)[cellID]),
    faceField_(cellID, mesh)
{
    if (!dofField_ || !mesh_)
    {
        FatalErrorInFunction
            << "Null pointer in GaussField(dofField,cellID,mesh)"
            << abort(FatalError);
    }

    // Build neighbor DOF list
    const labelList& ngh = mesh_->cells()[cellID_]->neighborCells();
    List<const cellDof<Type>*> cellsDof(ngh.size() + 1);

    cellsDof[0] = &(*dofField_)[cellID_];

    forAll(ngh, i)
    {
        const label nid = ngh[i];
        cellsDof[i + 1] = (nid >= 0 ? &(*dofField_)[nid] : nullptr);
    }

    faceField_.setCellsDof(cellsDof);
}


template<class Type>
GaussField<Type>::GaussField
(
    label cellID,
    const dgGeomMesh* mesh
)
:
    dofField_(nullptr),
    cellID_(cellID),
    mesh_(mesh),
    cellField_(cellID, mesh),
    faceField_(cellID, mesh)
{}


template<class Type>
GaussField<Type>::GaussField
(
    label cellID,
    const dgGeomMesh* mesh,
    const Type& initVal
)
:
    dofField_(nullptr),
    cellID_(cellID),
    mesh_(mesh),
    cellField_(cellID, mesh, initVal),
    faceField_(cellID, mesh, initVal)
{}


template<class Type>
GaussField<Type>::GaussField(const GaussField<Type>& other)
:
    dofField_(other.dofField_),
    cellID_(other.cellID_),
    mesh_(other.mesh_),
    cellField_(other.cellField_),
    faceField_(other.faceField_)
{}

/*
template<class Type>
GaussField<Type>::GaussField(const tmp<GaussField<Type>>& t)
{
    if (t.isTmp())
    {
        // Thực hiện Move assignment từ đối tượng bên trong t
        this->operator=(std::move(const_cast<GaussField<Type>&>(t())));
        t.clear();
    }
    else
    {
        *this = t(); // Copy nếu là biến cố định
    }
}

template<class Type>
GaussField<Type>::GaussField
(
    const tmp<GaussField<Type>>& t
)
:
    GaussField<Type>(t())
{
    t.clear();
}
*/

template<class Type>
GaussField<Type>::GaussField
(
    const tmp<GaussField<Type>>& t
)
{
    GaussField<Type>& tf =
        const_cast<tmp<GaussField<Type>>&>(t).ref();

    dofField_ = tf.dofField_;
    cellID_   = tf.cellID_;
    mesh_     = tf.mesh_;

    cellField_ = std::move(tf.cellField_);
    faceField_ = std::move(tf.faceField_);

    t.clear();
}

template<class Type>
GaussField<Type>::GaussField
(
    GaussField<Type>&& other
) noexcept
:
    dofField_(other.dofField_),
    cellID_(other.cellID_),
    mesh_(other.mesh_),
    cellField_(std::move(other.cellField_)),
    faceField_(std::move(other.faceField_))
{}

template<class Type>
tmp<GaussField<Type>> GaussField<Type>::New
(
    label cellID,
    const dgGeomMesh* mesh
)
{
    return tmp<GaussField<Type>>
    (
        new GaussField<Type>(cellID, mesh)
    );
}

// * * * * * * Assignment * * * * * //

template<class Type>
inline GaussField<Type>& GaussField<Type>::operator=
(
    const GaussField<Type>& other
)
{
    if (this != &other)
    {
        dofField_  = other.dofField_;
        cellID_    = other.cellID_;
        mesh_      = other.mesh_;
        cellField_ = other.cellField_;
        faceField_ = other.faceField_;
    }
    return *this;
}

template<class Type>
inline GaussField<Type>& GaussField<Type>::operator=
(
    const Type& value
)
{
    cellField_ = value;
    faceField_ = value;
    return *this;
}

template<class Type>
GaussField<Type>& Foam::GaussField<Type>::operator=
(
    GaussField<Type>&& other
) noexcept
{
    if (this != &other)
    {
        dofField_  = other.dofField_;
        cellID_    = other.cellID_;
        mesh_      = other.mesh_;

        cellField_ = std::move(other.cellField_);
        faceField_ = std::move(other.faceField_);
    }

    return *this;
}

template<class Type>
GaussField<Type>& Foam::GaussField<Type>::operator=
(
    const tmp<GaussField<Type>>& t
)
{
    /*
    *this = t();
    t.clear();
    return *this;
    */

    GaussField<Type>& tf =
        const_cast<tmp<GaussField<Type>>&>(t).ref();

    dofField_ = tf.dofField_;
    cellID_   = tf.cellID_;
    mesh_     = tf.mesh_;

    cellField_ = std::move(tf.cellField_);
    faceField_ = std::move(tf.faceField_);

    t.clear();

    return *this;
}

// * * * * * Access helpers * * * * * //

template<class Type>
inline const dgGeomMesh* GaussField<Type>::dgMesh() const
{
    if (!mesh_)
    {
        FatalErrorInFunction
            << "GaussField: mesh_ pointer is null"
            << abort(FatalError);
    }
    return mesh_;
}


template<class Type>
inline Type& GaussField<Type>::cellValueAt(label gpI)
{
    return cellField_[gpI];
}

template<class Type>
inline const Type& GaussField<Type>::cellValueAt(label gpI) const
{
    return cellField_[gpI];
}


template<class Type>
inline Type& GaussField<Type>::facePlusValueAt(label f, label gp)
{
    return faceField_.plusValueOnFace(f, gp);
}

template<class Type>
inline Type& GaussField<Type>::faceMinusValueAt(label f, label gp)
{
    return faceField_.minusValueOnFace(f, gp);
}

template<class Type>
inline const Type& GaussField<Type>::facePlusValueAt(label f, label gp) const
{
    return faceField_.plusValueOnFace(f, gp);
}

template<class Type>
inline const Type& GaussField<Type>::faceMinusValueAt(label f, label gp) const
{
    return faceField_.minusValueOnFace(f, gp);
}

// * * * * * * * * Member Functions * * * * * * * * //

//- Interpolate from DOF values
template<class Type>
void GaussField<Type>::interpolateFromDof()
{
    cellField_.interpolateFromDof();
    faceField_.interpolateFromDof();
}

} // End namespace Foam

#endif
/* ************************************************************************* */
