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

template<class Type>
cellGaussField<Type>::cellGaussField
(
    const label cellID,
    const dgGeomMesh& dgMesh,
    const List<Type>& values
)
:
    cellID_(cellID),
    dgMesh_(dgMesh),
    cell_(dgMesh_.cells()[cellID_]),
    nGauss_(cell_->gaussPoints().size()),
    nDof_(0),                 // no DOF provided
    dof_(0),                  // empty DOF list
    values_(values)           // copy input values
{
    if (values_.size() != nGauss_)
    {
        FatalErrorInFunction
            << "Size of input values (" << values_.size()
            << ") does not match number of Gauss points (" << nGauss_ << ")"
            << abort(FatalError);
    }
}


// Construct from DOF object
template<class Type>
cellGaussField<Type>::cellGaussField
(
    const cellDof<Type>& dof
)
:
    cellID_(dof.cellID()),
    dgMesh_(dof.mesh()),
    cell_(dof.cell()),
    nGauss_(cell_->gaussPoints().size()),
    nDof_(dof.nDof()),
    dof_(dof.dof()),
    values_(nGauss_) // allocate space for Gauss point values
{}


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


// Initial value constructor
template<class Type>
Foam::cellGaussField<Type>::cellGaussField
(
    const label cellID,
    const dgGeomMesh& dgMesh,
    const Type& initVal
)
:
    cellID_(cellID),
    dgMesh_(dgMesh),
    cell_(dgMesh_.cells()[cellID_]),
    nGauss_(cell_->gaussPoints().size()),
    nDof_(0),
    dof_(0),
    values_(nGauss_, initVal) // fill with initVal
{}

// Constructor without initial values
// This is useful when you want to create an empty field and fill it later
template<class Type>
Foam::cellGaussField<Type>::cellGaussField
(
    const label cellID,
    const dgGeomMesh& dgMesh
)
:
    cellID_(cellID),
    dgMesh_(dgMesh),
    cell_(dgMesh_.cells()[cellID_]),
    nGauss_(cell_->gaussPoints().size()),
    nDof_(0),
    dof_(0),
    values_(nGauss_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void cellGaussField<Type>::print() const
{
    Info << "Gauss field values for cell " << cellID_ << nl;

    for (label i = 0; i < nGauss_; ++i)
    {
        Info << "  Gauss pt " << i << ": " << values_[i] << nl;
    }
}

template<class Type>
void Foam::cellGaussField<Type>::interpolateFromDof()
{
    // Access Gauss points (not used here but can be helpful later)
    const List<vector>& gaussPts = cell_->gaussPoints();

    // Access basis function values at Gauss points
    const List<List<scalar>>& b = cell_->basis();

    // Loop over each Gauss point
    for (label gp = 0; gp < nGauss_; ++gp)
    {
        // Initialize interpolated value
        Type val = pTraits<Type>::zero; // zero for scalar, vector, tensor... this is OpenFOAM's standard

        // Sum over each DOF
        for (label k = 0; k < nDof_; ++k)
        {
            val += dof_[k] * b[gp][k];
        }

        // Store interpolated value
        values_[gp] = val;
    }
}

template<class Type>
Foam::cellGaussField<Type>& Foam::cellGaussField<Type>::operator=
(
    const cellGaussField<Type>& other
)
{
    // Self-assignment check
    if (this == &other)
    {
        return *this;
    }

    // Check if same cell
    if (cellID_ != other.cellID_)
    {
        FatalErrorInFunction
            << "Assignment between different cells is not allowed. "
            << "cellID_ = " << cellID_ << ", other.cellID_ = " << other.cellID_
            << abort(FatalError);
    }

    // Copy values (same size assumed)
    values_ = other.values_;
    dof_    = other.dof_;
    nDof_   = other.nDof_;

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

} // End namespace Foam

// ************************************************************************* //
