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
    Foam::faceGaussField

\*---------------------------------------------------------------------------*/

// faceGaussFieldVectorI.H
// Implementation of template specializations for faceGaussField<vector>

namespace Foam
{

//- Limit all Gauss-point values (plus/minus) to a maximum value
template<>
inline void faceGaussField<scalar>::limitMax(const scalar maxVal)
{
    for (label gp = 0; gp < nGauss(); ++gp)
    {
        if (plusValues_[gp] > maxVal)
        {
            plusValues_[gp] = maxVal;
        }
        if (minusValues_[gp] > maxVal)
        {
            minusValues_[gp] = maxVal;
        }
    }
}


//- Limit all Gauss-point values (plus/minus) to a minimum value
template<>
inline void faceGaussField<scalar>::limitMin(const scalar minVal)
{
    for (label gp = 0; gp < nGauss(); ++gp)
    {
        if (plusValues_[gp] < minVal)
        {
            plusValues_[gp] = minVal;
        }
        if (minusValues_[gp] < minVal)
        {
            minusValues_[gp] = minVal;
        }
    }
}


//- Limit all Gauss-point values (plus/minus) to lie within [minVal, maxVal]
template<>
inline void faceGaussField<scalar>::limitRange
(
    const scalar minVal,
    const scalar maxVal
)
{
    for (label gp = 0; gp < nGauss(); ++gp)
    {
        if (plusValues_[gp] < minVal)
        {
            plusValues_[gp] = minVal;
        }
        else if (plusValues_[gp] > maxVal)
        {
            plusValues_[gp] = maxVal;
        }

        if (minusValues_[gp] < minVal)
        {
            minusValues_[gp] = minVal;
        }
        else if (minusValues_[gp] > maxVal)
        {
            minusValues_[gp] = maxVal;
        }
    }
}

} // namespace Foam

// ************************************************************************* //

