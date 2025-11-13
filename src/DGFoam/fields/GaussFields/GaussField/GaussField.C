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

#include "GaussField.H"
#include "IOstreams.H"
namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::GaussField<Type>::GaussField()
:
    dofField_(nullptr),
    cellID_(-1),
    mesh_(nullptr),
    cellField_(),
    faceField_()
{}

// Construct from dofField, cellID and mesh
template<class Type>
Foam::GaussField<Type>::GaussField
(
    const dofField<Type>* dofField,
    label cellID,
    const dgGeomMesh* mesh
)
:
    dofField_(dofField),
    cellID_(cellID),
    mesh_(mesh),
    cellField_(mesh, &(*dofField_)[cellID_]),
    faceField_(cellID,mesh)
{
    if (!dofField_ || !mesh_)
    {
        FatalErrorInFunction
            << "Null dofField_ or mesh_ pointer in GaussField constructor"
            << abort(FatalError);
    }

    const labelList& neighbors = mesh_->cells()[cellID_]->neighborCells();
    List<const cellDof<Type>*> cellsDof(neighbors.size()+1);

    // The first entry is the cell itself
    cellsDof[0] = &(*dofField_)[cellID_];

    forAll(neighbors, i)
    {
        const label nid = neighbors[i];
        if (nid >= 0)
        {
            cellsDof[i+1] = &(*dofField_)[nid];
        }
        else
        {
            cellsDof[i+1] = nullptr;
        }
    }
    faceField_.setCellsDof(cellsDof);
}


// Construct empty container only

template<class Type>
Foam::GaussField<Type>::GaussField
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


// Construct container with initial value

template<class Type>
Foam::GaussField<Type>::GaussField
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

// Copy constructor
template<class Type>
Foam::GaussField<Type>::GaussField
(
    const GaussField<Type>& other
)
:
    dofField_(other.dofField_),
    cellID_(other.cellID_),
    mesh_(other.mesh_),
    cellField_(other.cellField_),
    faceField_(other.faceField_)
{}

// * * * * * * * * * * * * * * * Assignment * * * * * * * * * * * * * * * * //

template<class Type>
Foam::GaussField<Type>& Foam::GaussField<Type>::operator=
(
    const GaussField<Type>& other
)
{
    if (this == &other)
    {
        return *this;
    }

    dofField_ = other.dofField_;
    cellID_ = other.cellID_;
    cellField_ = other.cellField_;
    faceField_ = other.faceField_;
    mesh_ = other.mesh_;

    return *this;
}

template<class Type>
Foam::GaussField<Type>& Foam::GaussField<Type>::operator=
(
    const Type& value
)
{
    cellField_ = value;
    faceField_ = value;

    return *this;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::GaussField<Type>::interpolateFromDof()
{
    cellField_.interpolateFromDof();
    faceField_.interpolateFromDof();
}

// * * * * * * * * * * * * * * Template Instantiations  * * * * * * * * * * //

template class Foam::GaussField<Foam::scalar>;
template class Foam::GaussField<Foam::vector>;
template class Foam::GaussField<Foam::tensor>;
template class Foam::GaussField<Foam::symmTensor>;
template class Foam::GaussField<Foam::sphericalTensor>;

// ************************************************************************* //
} // End namespace Foam
