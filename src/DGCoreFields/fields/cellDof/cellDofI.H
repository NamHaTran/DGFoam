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
    Foam::cellDof<Type>

Description
    Inline implementations for the Foam::cellDof class.

SourceFiles
    cellDofI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellDofI_H
#define Foam_cellDofI_H

namespace Foam
{

// * * * * * * * * * * * Constructors  * * * * * * * * * * * * //

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
    const tmp<cellDof<Type>>& t
)
:
    cellID_(-1),
    nDof_(0),
    dof_()
{
    cellDof<Type>& tf = const_cast<tmp<cellDof<Type>>&>(t).ref();

    cellID_ = tf.cellID_;
    nDof_   = tf.nDof_;
    dof_.transfer(tf.dof_);

    t.clear();
}


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


// * * * * * * * * * * * Assignment  * * * * * * * * * * * * * //

template<class Type>
cellDof<Type>& cellDof<Type>::operator=(const cellDof<Type>& other)
{
    if (this != &other)
    {
        cellID_ = other.cellID_;
        nDof_ = other.nDof_;
        dof_ = other.dof_;
    }

    return *this;
}


template<class Type>
cellDof<Type>& cellDof<Type>::operator=
(
    const tmp<cellDof<Type>>& t
)
{
    cellDof<Type>& tf = const_cast<tmp<cellDof<Type>>&>(t).ref();

    cellID_ = tf.cellID_;
    nDof_   = tf.nDof_;
    dof_.transfer(tf.dof_);

    t.clear();

    return *this;
}

} // End namespace Foam

#endif

// ************************************************************************* //
