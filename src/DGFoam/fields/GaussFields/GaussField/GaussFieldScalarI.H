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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Limit all Gauss-point values (cell and face) to a maximum value
template<>
inline void GaussField<scalar>::limitMax(const scalar maxVal)
{
    cellField_.limitMax(maxVal);
    faceField_.limitMax(maxVal);
}


//- Limit all Gauss-point values (cell and face) to a minimum value
template<>
inline void GaussField<scalar>::limitMin(const scalar minVal)
{
    cellField_.limitMin(minVal);
    faceField_.limitMin(minVal);
}


//- Limit all Gauss-point values (cell and face) to lie within [minVal, maxVal]
template<>
inline void GaussField<scalar>::limitRange
(
    const scalar minVal,
    const scalar maxVal
)
{
    // Delegate directly to subfields
    cellField_.limitRange(minVal, maxVal);
    faceField_.limitRange(minVal, maxVal);
}

} // End namespace Foam

// ************************************************************************* //

