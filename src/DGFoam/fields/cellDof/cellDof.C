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

#include "cellDof.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cellDof<Type>::cellDof()
:
    cellID_(-1),
    nDof_(0),
    dof_()
{}


template<class Type>
cellDof<Type>::cellDof(const cellDof<Type>& other)
:
    cellID_(other.cellID_),
    nDof_(other.nDof_),
    dof_(other.dof_)
{}


template<class Type>
cellDof<Type>::cellDof
(
    const label cellID,
    const label nDof,
    const UList<Type>& inputDof
)
:
    cellID_(cellID),
    nDof_(nDof),
    dof_(nDof, pTraits<Type>::zero)
{
    const label inputSize = inputDof.size();
    const label minSize = min(inputSize, nDof_);

    for (label i = 0; i < minSize; ++i)
    {
        dof_[i] = inputDof[i];
    }
}


// * * * * * * * * * * * * * * * * Assignment * * * * * * * * * * * * * * * * //

template<class Type>
cellDof<Type>& cellDof<Type>::operator=(const cellDof<Type>& other)
{
    if (this != &other)
    {
        cellID_ = other.cellID_;
        nDof_   = other.nDof_;
        dof_    = other.dof_;
    }

    return *this;
}


// * * * * * * * * * * * * * * * * Instantiation * * * * * * * * * * * * * * //

template class cellDof<scalar>;
template class cellDof<vector>;
template class cellDof<tensor>;
template class cellDof<symmTensor>;
template class cellDof<sphericalTensor>;

} // End namespace Foam

// ************************************************************************* //
