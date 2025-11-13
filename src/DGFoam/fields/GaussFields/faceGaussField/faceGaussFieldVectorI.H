/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR by YOUR NAME
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

template<>
inline faceGaussField<scalar> faceGaussField<vector>::x() const
{
    faceGaussField<scalar> result(cellID_, mesh_);
    forAll(plusValues_, i)
    {
        result.plusValueAt(i)  = plusValues_[i].x();
        result.minusValueAt(i) = minusValues_[i].x();
    }
    return result;
}

template<>
inline faceGaussField<scalar> faceGaussField<vector>::y() const
{
    faceGaussField<scalar> result(cellID_, mesh_);
    forAll(plusValues_, i)
    {
        result.plusValueAt(i)  = plusValues_[i].y();
        result.minusValueAt(i) = minusValues_[i].y();
    }
    return result;
}

template<>
inline faceGaussField<scalar> faceGaussField<vector>::z() const
{
    faceGaussField<scalar> result(cellID_, mesh_);
    forAll(plusValues_, i)
    {
        result.plusValueAt(i)  = plusValues_[i].z();
        result.minusValueAt(i) = minusValues_[i].z();
    }
    return result;
}

template<>
inline void faceGaussField<vector>::setX(const faceGaussField<scalar>& fx)
{
    /*
    if (nGauss_ != fx.nGauss_ || cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Incompatible field dimensions or cell mismatch in setX()"
            << abort(FatalError);
    }
    */
    forAll(plusValues_, i)
    {
        plusValues_[i].x()  = fx.plusValue(i);
        minusValues_[i].x() = fx.minusValue(i);
    }
}

template<>
inline void faceGaussField<vector>::setY(const faceGaussField<scalar>& fy)
{
    /*
    if (nGauss_ != fx.nGauss_ || cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Incompatible field dimensions or cell mismatch in setX()"
            << abort(FatalError);
    }
    */
    forAll(plusValues_, i)
    {
        plusValues_[i].y()  = fy.plusValue(i);
        minusValues_[i].y() = fy.minusValue(i);
    }
}

template<>
inline void faceGaussField<vector>::setZ(const faceGaussField<scalar>& fz)
{
    /*
    if (nGauss_ != fx.nGauss_ || cellID_ != fx.cellID_)
    {
        FatalErrorInFunction
            << "Incompatible field dimensions or cell mismatch in setX()"
            << abort(FatalError);
    }
    */
    forAll(plusValues_, i)
    {
        plusValues_[i].z()  = fz.plusValue(i);
        minusValues_[i].z() = fz.minusValue(i);
    }
}

template<>
inline void faceGaussField<vector>::setXYZ
(
    const faceGaussField<scalar>& fx,
    const faceGaussField<scalar>& fy,
    const faceGaussField<scalar>& fz
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

    forAll(plusValues_, i)
    {
        plusValues_[i]  = vector(fx.plusValue(i),  fy.plusValue(i),  fz.plusValue(i));
        minusValues_[i] = vector(fx.minusValue(i), fy.minusValue(i), fz.minusValue(i));
    }
}

//- Limit vector magnitudes (plus/minus) to a maximum value
template<>
inline void faceGaussField<vector>::limitMax(const scalar magMax)
{
    for (label gp = 0; gp < nGauss(); ++gp)
    {
        vector& vPlus = plusValues_[gp];
        vector& vMinus = minusValues_[gp];

        const scalar magP = mag(vPlus);
        const scalar magM = mag(vMinus);

        if (magP > magMax && magP > VSMALL)
        {
            vPlus *= (magMax / magP);
        }

        if (magM > magMax && magM > VSMALL)
        {
            vMinus *= (magMax / magM);
        }
    }
}


//- Limit vector magnitudes (plus/minus) to a minimum value
template<>
inline void faceGaussField<vector>::limitMin(const scalar magMin)
{
    for (label gp = 0; gp < nGauss(); ++gp)
    {
        vector& vPlus = plusValues_[gp];
        vector& vMinus = minusValues_[gp];

        const scalar magP = mag(vPlus);
        const scalar magM = mag(vMinus);

        if (magP < magMin && magP > VSMALL)
        {
            vPlus *= (magMin / magP);
        }

        if (magM < magMin && magM > VSMALL)
        {
            vMinus *= (magMin / magM);
        }
    }
}


//- Limit vector magnitudes (plus/minus) to lie within [magMin, magMax]
template<>
inline void faceGaussField<vector>::limitRange
(
    const scalar magMin,
    const scalar magMax
)
{
    for (label gp = 0; gp < nGauss(); ++gp)
    {
        vector& vPlus = plusValues_[gp];
        vector& vMinus = minusValues_[gp];

        const scalar magP = mag(vPlus);
        const scalar magM = mag(vMinus);

        // Limit plus side
        if (magP > magMax && magP > VSMALL)
        {
            vPlus *= (magMax / magP);
        }
        else if (magP < magMin && magP > VSMALL)
        {
            vPlus *= (magMin / magP);
        }

        // Limit minus side
        if (magM > magMax && magM > VSMALL)
        {
            vMinus *= (magMax / magM);
        }
        else if (magM < magMin && magM > VSMALL)
        {
            vMinus *= (magMin / magM);
        }
    }
}

} // namespace Foam

// ************************************************************************* //

