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

#include "cellGaussField.H"
#include "vector.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Default constructor
template<class Type>
Foam::cellGaussField<Type>::cellGaussField()
:
    cellID_(-1),
    dgMesh_(nullptr),
    cell_(nullptr),
    nGauss_(0),
    nDof_(0),
    dof_(nullptr),
    values_()
{}

// Construct from DOF object
template<class Type>
cellGaussField<Type>::cellGaussField
(
    const dgGeomMesh* dgMesh,
    const cellDof<Type>* dof
)
:
    cellID_(dof->cellID()),
    dgMesh_(dgMesh),
    cell_(dgMesh->cells()[cellID_]),
    nGauss_(cell_->gaussPoints().size()),
    nDof_(dof->nDof()),
    dof_(dof),
    values_(nGauss_)
{
    if (!dgMesh)
    {
        FatalErrorInFunction << "dgMesh pointer is null" << abort(FatalError);
    }
}


template<class Type>
Foam::cellGaussField<Type>::cellGaussField
(
    const cellGaussField<Type>& other
)
:
    cellID_(other.cellID_),       // Copy cell ID
    dgMesh_(other.dgMesh_),       // Copy pointer to mesh (non-owning)
    cell_(other.cell_),           // Copy pointer to cell (non-owning)
    nGauss_(other.nGauss_),       // Copy number of Gauss points
    nDof_(other.nDof_),           // Copy number of DOFs
    dof_(other.dof_),             // Copy pointer to DOF (non-owning)
    values_(other.values_)        // Deep copy of Gauss-point values
{
    // Nothing else to do
    // All pointers are non-owning, so shallow copy is safe
}


// Initial value constructor
template<class Type>
Foam::cellGaussField<Type>::cellGaussField
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
    if (!dgMesh)
    {
        FatalErrorInFunction << "dgMesh pointer is null" << abort(FatalError);
    }
}


// Constructor without initial values
template<class Type>
Foam::cellGaussField<Type>::cellGaussField
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
    if (!dgMesh)
    {
        FatalErrorInFunction << "dgMesh pointer is null" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::cellGaussField<Type>::interpolateFromDof()
{
    //const List<vector>& gaussPts = cell_.gaussPoints();
    const List<List<scalar>>& b = cell_->basis();

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        Type val = pTraits<Type>::zero;

        for (label k = 0; k < nDof_; ++k)
        {
            val += dof_->dof()[k] * b[gp][k];
        }

        values_[gp] = val;
    }
}

template<class Type>
Foam::cellGaussField<Type>& Foam::cellGaussField<Type>::operator=
(
    const cellGaussField<Type>& other
)
{
    if (this == &other)
    {
        return *this;
    }

    // Copy basic attributes
    cellID_ = other.cellID_;           // Copy cell ID
    nGauss_ = other.nGauss_;           // Copy number of Gauss points
    nDof_   = other.nDof_;             // Copy number of DOFs

    // Copy pointer references (non-owning pointers)
    dgMesh_ = other.dgMesh_;           // Copy pointer to mesh (no ownership)
    cell_   = other.cell_;             // Copy pointer to cell (no ownership)
    dof_    = other.dof_;              // Copy pointer to DOF (no ownership)

    // Deep copy of values (field data at Gauss points)
    values_ = other.values_;

    return *this;
}

template<class Type>
Foam::cellGaussField<Type>& Foam::cellGaussField<Type>::operator=
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

// * * * * * * * * * * * * * * Template instantiation * * * * * * * * * * * * //
template class Foam::cellGaussField<scalar>;
template class Foam::cellGaussField<vector>;
template class Foam::cellGaussField<tensor>;
template class Foam::cellGaussField<symmTensor>;
template class Foam::cellGaussField<sphericalTensor>;

} // End namespace Foam

// ************************************************************************* //

