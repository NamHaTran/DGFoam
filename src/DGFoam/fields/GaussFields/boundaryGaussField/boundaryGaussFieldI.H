/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2026 Ha Nam Tran
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

#ifndef Foam_boundaryGaussFieldI_H
#define Foam_boundaryGaussFieldI_H

#include "boundaryGaussField.H"

namespace Foam
{

template<class Type>
inline boundaryGaussField<Type>::boundaryGaussField()
:
    values_()
{}


template<class Type>
inline boundaryGaussField<Type>::boundaryGaussField
(
    const label nSize
)
:
    values_(nSize)
{}


template<class Type>
inline boundaryGaussField<Type>::boundaryGaussField
(
    const boundaryGaussField<Type>& other
)
:
    values_(other.values_)
{}


template<class Type>
inline boundaryGaussField<Type>& boundaryGaussField<Type>::operator=
(
    const boundaryGaussField<Type>& other
)
{
    if (this != &other)
    {
        values_ = other.values_;
    }

    return *this;
}


template<class Type>
inline boundaryGaussField<Type>& boundaryGaussField<Type>::operator=
(
    const Type& value
)
{
    forAll(values_, i)
    {
        values_[i] = value;
    }

    return *this;
}


template<class Type>
inline label boundaryGaussField<Type>::size() const
{
    return values_.size();
}


template<class Type>
inline Type& boundaryGaussField<Type>::operator[](const label i)
{
    return values_[i];
}


template<class Type>
inline const Type& boundaryGaussField<Type>::operator[](const label i) const
{
    return values_[i];
}

} // End namespace Foam

#endif

// ************************************************************************* //
