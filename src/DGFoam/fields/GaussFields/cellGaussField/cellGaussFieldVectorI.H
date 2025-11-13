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

Class
    Foam::cellGaussField

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
inline cellGaussField<scalar> cellGaussField<vector>::x() const
{
    cellGaussField<scalar> result(cellID_, dgMesh_);
    forAll(values_, i)
    {
        result[i] = values_[i].x();
    }
    return result;
}

template<>
inline cellGaussField<scalar> cellGaussField<vector>::y() const
{
    cellGaussField<scalar> result(cellID_, dgMesh_);
    forAll(values_, i)
    {
        result[i] = values_[i].y();
    }
    return result;
}

template<>
inline cellGaussField<scalar> cellGaussField<vector>::z() const
{
    cellGaussField<scalar> result(cellID_, dgMesh_);
    forAll(values_, i)
    {
        result[i] = values_[i].z();
    }
    return result;
}

//- Set x-component from scalar field
template<>
inline void cellGaussField<vector>::setX(const cellGaussField<scalar>& fx)
{
    /*
    if (nGauss_ != fx.nGauss_ || cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Incompatible field dimensions or cell mismatch in setX()"
            << abort(FatalError);
    }
    */

    forAll(values_, i)
    {
        values_[i].x() = fx[i];
    }
}

//- Set y-component from scalar field
template<>
inline void cellGaussField<vector>::setY(const cellGaussField<scalar>& fy)
{
    /*
    if (nGauss_ != fx.nGauss_ || cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Incompatible field dimensions or cell mismatch in setX()"
            << abort(FatalError);
    }
    */

    forAll(values_, i)
    {
        values_[i].y() = fy[i];
    }
}

//- Set z-component from scalar field
template<>
inline void cellGaussField<vector>::setZ(const cellGaussField<scalar>& fz)
{
    /*
    if (nGauss_ != fx.nGauss_ || cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Incompatible field dimensions or cell mismatch in setX()"
            << abort(FatalError);
    }
    */

    forAll(values_, i)
    {
        values_[i].z() = fz[i];
    }
}


//- Set all components from 3 scalar fields
template<>
inline void cellGaussField<vector>::setXYZ
(
    const cellGaussField<scalar>& fx,
    const cellGaussField<scalar>& fy,
    const cellGaussField<scalar>& fz
)
{
    /*
    if
    (
        nGauss_ != fx.nGauss_ ||
        nGauss_ != fy.nGauss_ ||
        nGauss_ != fz.nGauss_ ||
        cellID_ != fx.cellID_ ||
        cellID_ != fy.cellID_ ||
        cellID_ != fz.cellID_
    )
    {
        FatalErrorInFunction
            << "Incompatible field dimensions or cell mismatch in setXYZ()"
            << abort(FatalError);
    }
    */

    forAll(values_, i)
    {
        values_[i] = vector(fx[i], fy[i], fz[i]);
    }
}

//- Limit vector magnitudes to a maximum value
template<>
inline void cellGaussField<vector>::limitMax(const scalar magMax)
{
    for (label gp = 0; gp < size(); ++gp)
    {
        vector& v = values_[gp];
        const scalar m = mag(v);

        if (m > magMax && m > VSMALL)
        {
            // Scale vector down to magnitude = magMax
            v *= (magMax / m);
        }
    }
}


//- Limit vector magnitudes to a minimum value
template<>
inline void cellGaussField<vector>::limitMin(const scalar magMin)
{
    for (label gp = 0; gp < size(); ++gp)
    {
        vector& v = values_[gp];
        const scalar m = mag(v);

        if (m < magMin && m > VSMALL)
        {
            // Scale vector up to magnitude = magMin
            v *= (magMin / m);
        }
    }
}


//- Limit vector magnitudes to lie within [magMin, magMax]
template<>
inline void cellGaussField<vector>::limitRange
(
    const scalar magMin,
    const scalar magMax
)
{
    for (label gp = 0; gp < size(); ++gp)
    {
        vector& v = values_[gp];
        const scalar m = mag(v);

        if (m > magMax && m > VSMALL)
        {
            v *= (magMax / m);
        }
        else if (m < magMin && m > VSMALL)
        {
            v *= (magMin / m);
        }
    }
}

// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //

