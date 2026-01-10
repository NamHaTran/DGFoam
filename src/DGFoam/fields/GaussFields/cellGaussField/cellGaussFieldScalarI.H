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

Class
    Foam::cellGaussField

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Limit all Gauss-point values to a maximum value
template<>
inline void cellGaussField<scalar>::limitMax(const scalar Tmax)
{
    for (label gp = 0; gp < size(); ++gp)
    {
        if (values_[gp] > Tmax)
        {
            values_[gp] = Tmax;
        }
    }
}


//- Limit all Gauss-point values to a minimum value
template<>
inline void cellGaussField<scalar>::limitMin(const scalar Tmin)
{
    for (label gp = 0; gp < size(); ++gp)
    {
        if (values_[gp] < Tmin)
        {
            values_[gp] = Tmin;
        }
    }
}


//- Limit all Gauss-point values to lie within [Tmin, Tmax]
template<>
inline void cellGaussField<scalar>::limitRange
(
    const scalar Tmin,
    const scalar Tmax
)
{
    for (label gp = 0; gp < size(); ++gp)
    {
        if (values_[gp] < Tmin)
        {
            values_[gp] = Tmin;
        }
        else if (values_[gp] > Tmax)
        {
            values_[gp] = Tmax;
        }
    }
}

// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //

