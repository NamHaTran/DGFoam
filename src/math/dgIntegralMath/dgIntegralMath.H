/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Description
    Mathematical utilities for calculating integral.
SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgIntegralMath_H
#define Foam_dgIntegralMath_H

#include "cellGaussField.H"
#include "faceGaussField.H"
#include "dgGeomFace.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Volume integral:  ∫ F dV  ≈ Σ_g F(g) w(g)
\*---------------------------------------------------------------------------*/
inline scalar volumeIntegral(const cellGaussField<vector>& gradPhi, const cellGaussField<vector>& F)
{
    scalar volInt = Zero;

    const label nGauss = F.size();
    const List<scalar>& w = F.cell()->weights();

    // Cell has straight edges/faces, so Jacobian is constant (= cell volume)
    const scalar cellVol = F.cell()->volume();

    for (label g = 0; g < nGauss; ++g)
    {
        volInt += cellVol * (gradPhi[g] & F[g]) * w[g];
    }

    return volInt;
}

inline vector volumeIntegral(const cellGaussField<vector>& gradPhi, const cellGaussField<tensor>& F)
{
    vector volInt = Zero;

    const label nGauss = F.size();
    const List<scalar>& w = F.cell()->weights();

    // Cell has straight edges/faces, so Jacobian is constant (= cell volume)
    const scalar cellVol = F.cell()->volume();

    for (label g = 0; g < nGauss; ++g)
    {
        volInt += cellVol * (gradPhi[g] & F[g]) * w[g];
    }

    return volInt;
}

/*---------------------------------------------------------------------------*\
     Surface integral for vector flux:
       F is vector, surface integral returns scalar
       Integral = ∫ phi (F · n) dS ≈ Σ phi (F(g) & n) w(g)
\*---------------------------------------------------------------------------*/
inline scalar surfaceIntegral(const faceGaussField<scalar>& phi, const faceGaussField<vector>& F)
{
    scalar surfInt = 0.0;

    const label nFaces = F.nFaces();

    for (label faceI = 0; faceI < nFaces; ++faceI)
    {
        const vector n = F.faces()[faceI]->normal();
        const label nGauss = F.nGaussPerFace();
        const List<scalar>& w = F.faces()[faceI]->weights();

        // Face has straight edges, so Jacobian is constant (= face area)
        const scalar faceArea = F.faces()[faceI]->area();

        for (label g = 0; g < nGauss; ++g)
        {
            // F(g) is vector ⇒ dot with n to get scalar flux
            surfInt += faceArea * phi.minusValueOnFace(faceI, g) * (F.fluxOnFace(faceI, g) & n) * w[g];
        }
    }

    return surfInt;
}


/*---------------------------------------------------------------------------*\
     Surface integral for tensor flux:
       F is tensor, but F.fluxOnFace already returns T·n (a vector)
       Surface integral returns vector:
         ∫ (T·n) dS ≈ Σ (T·n)(g) w(g)
\*---------------------------------------------------------------------------*/
inline vector surfaceIntegral(const faceGaussField<scalar>& phi, const faceGaussField<tensor>& F)
{
    vector surfInt = Zero;

    const label nFaces = F.nFaces();

    for (label faceI = 0; faceI < nFaces; ++faceI)
    {
        const label nGauss = F.nGaussPerFace();
        const List<scalar>& w = F.faces()[faceI]->weights();

        // Face has straight edges, so Jacobian is constant (= face area)
        const scalar faceArea = F.faces()[faceI]->area();

        for (label g = 0; g < nGauss; ++g)
        {
            // fluxOnFace already returns T·n (vector)
            surfInt += faceArea * phi.minusValueOnFace(faceI, g) * F.fluxOnFace(faceI, g) * w[g];
        }
    }

    return surfInt;
}

} // End namespace Foam

#endif

// ************************************************************************* //

