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

Class
    Foam::faceGaussField

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return X-component GaussField<scalar>
template<>
inline GaussField<scalar> GaussField<vector>::x() const
{
    GaussField<scalar> result(cellID_, mesh_);

    result.cellField() = cellField_.x();
    result.faceField() = faceField_.x();

    return result;
}


//- Return Y-component GaussField<scalar>
template<>
inline GaussField<scalar> GaussField<vector>::y() const
{
    GaussField<scalar> result(cellID_, mesh_);

    result.cellField() = cellField_.y();
    result.faceField() = faceField_.y();

    return result;
}


//- Return Z-component GaussField<scalar>
template<>
inline GaussField<scalar> GaussField<vector>::z() const
{
    GaussField<scalar> result(cellID_, mesh_);

    result.cellField() = cellField_.z();
    result.faceField() = faceField_.z();

    return result;
}

//- Set X-component from GaussField<scalar>
template<>
inline void GaussField<vector>::setX(const GaussField<scalar>& fx)
{
    /*
    if (cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Cell ID mismatch in setX()" << abort(FatalError);
    }
    */

    cellField_.setX(fx.cellField());
    faceField_.setX(fx.faceField());
}


//- Set Y-component from GaussField<scalar>
template<>
inline void GaussField<vector>::setY(const GaussField<scalar>& fy)
{
    /*
    if (cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Cell ID mismatch in setX()" << abort(FatalError);
    }
    */

    cellField_.setY(fy.cellField());
    faceField_.setY(fy.faceField());
}


//- Set Z-component from GaussField<scalar>
template<>
inline void GaussField<vector>::setZ(const GaussField<scalar>& fz)
{
    /*
    if (cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Cell ID mismatch in setX()" << abort(FatalError);
    }
    */

    cellField_.setZ(fz.cellField());
    faceField_.setZ(fz.faceField());
}


//- Set all components from 3 scalar GaussFields
template<>
inline void GaussField<vector>::setXYZ
(
    const GaussField<scalar>& fx,
    const GaussField<scalar>& fy,
    const GaussField<scalar>& fz
)
{
    /*
    if (cellID_ != fx.cellID_ || cellID_ != fy.cellID_ || cellID_ != fz.cellID_)
    {
        FatalErrorInFunction
            << "Cell ID mismatch in setXYZ()" << abort(FatalError);
    }
    */

    cellField_.setXYZ(fx.cellField(), fy.cellField(), fz.cellField());
    faceField_.setXYZ(fx.faceField(), fy.faceField(), fz.faceField());
}

//- Limit all Gauss-point values (cell and face) to a maximum value
template<>
inline void GaussField<vector>::limitMax(const scalar maxVal)
{
    cellField_.limitMax(maxVal);
    faceField_.limitMax(maxVal);
}


//- Limit all Gauss-point values (cell and face) to a minimum value
template<>
inline void GaussField<vector>::limitMin(const scalar minVal)
{
    cellField_.limitMin(minVal);
    faceField_.limitMin(minVal);
}


//- Limit all Gauss-point values (cell and face) to lie within [minVal, maxVal]
template<>
inline void GaussField<vector>::limitRange
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

