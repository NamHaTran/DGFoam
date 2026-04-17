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

Description
    Mathematical utilities for calculating integral.
SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgIntegralMath_H
#define Foam_dgIntegralMath_H

#include "cellGaussField.H"
#include "faceGaussField.H"
#include "dgGeomFace.H"
#include "products.H"

namespace Foam
{

/**
 * \brief Volume integral for vector field F, used in DG weak form. F is vector
 * 
 * \param gradPhi Gradient of basis function at Gauss points
 * \param F       Vector field at Gauss points
 */
inline scalar volumeIntegral(const cellGaussField<vector>& gradPhi, const cellGaussField<vector>& F)
{
    scalar volInt = Zero;

    const label nGauss = F.size();
    const List<scalar>& w = F.cell()->weights();
    const List<scalar>& J3D = F.cell()->J3D();

    for (label g = 0; g < nGauss; ++g)
    {
        volInt += J3D[g] * (gradPhi[g] & F[g]) * w[g];
    }

    return volInt;
}

/**
 * \brief Volume integral for vector field F, used in DG weak form. F is tensor
 * 
 * \param gradPhi Gradient of basis function at Gauss points
 * \param F       Vector field at Gauss points
 */
inline vector volumeIntegral(const cellGaussField<vector>& gradPhi, const cellGaussField<tensor>& F)
{
    vector volInt = Zero;

    const label nGauss = F.size();
    const List<scalar>& w = F.cell()->weights();
    const List<scalar>& J3D = F.cell()->J3D();

    for (label g = 0; g < nGauss; ++g)
    {
        volInt += J3D[g] * (gradPhi[g] & F[g]) * w[g];
    }

    return volInt;
}

/**
 * \brief Volume integral for the gradient-projection weak form.
 *
 * Computes \f$\int_\Omega C\,\nabla\phi\,dV\f$.
 * If C is scalar the result is vector; if C is vector the result is tensor.
 *
 * \param gradPhi Gradient of basis function at Gauss points
 * \param C       Projected field at cell Gauss points
 */
template<class Type>
inline typename outerProduct<Type, vector>::type gradProjectionVolumeIntegral
(
    const cellGaussField<vector>& gradPhi,
    const cellGaussField<Type>& C
)
{
    typename outerProduct<Type, vector>::type volInt =
        pTraits<typename outerProduct<Type, vector>::type>::zero;

    const label nGauss = C.size();
    const List<scalar>& w = C.cell()->weights();
    const List<scalar>& J3D = C.cell()->J3D();

    for (label g = 0; g < nGauss; ++g)
    {
        volInt += J3D[g] * (C[g] * gradPhi[g]) * w[g];
    }

    return volInt;
}

/**
 * \brief General volume integral for scalar  field F over the cell
 * 
 * \param F Scalar field at Gauss points
 * \param J Jacobian determinants at Gauss points
 * \param w Weights at Gauss points
 */
inline scalar volumeIntegral(const List<scalar>& F, const List<scalar>& J, const List<scalar>& w)
{
    scalar volInt = 0;

    const label nGauss = F.size();

    for (label g = 0; g < nGauss; ++g)
    {
        volInt += J[g] * F[g] * w[g];
    }

    return volInt;
}

/**
 * \brief General volume integral for vector field F over the cell
 * 
 * \param F Vector field at Gauss points
 * \param J Jacobian determinants at Gauss points
 * \param w Weights at Gauss points
 */
inline vector volumeIntegral(const List<vector>& F, const List<scalar>& J, const List<scalar>& w)
{
    vector volInt = Zero;

    const label nGauss = F.size();

    for (label g = 0; g < nGauss; ++g)
    {
        volInt += J[g] * F[g] * w[g];
    }

    return volInt;
}

/**
 * \brief Surface integral for vector flux
 * 
 * \param phi Scalar field at face Gauss points
 * \param F   Vector field at face Gauss points
 */
inline scalar surfaceIntegral(const faceGaussField<scalar>& phi, const faceGaussField<vector>& F)
{
    scalar surfInt = 0.0;

    const label nFaces = F.nFaces();

    for (label faceI = 0; faceI < nFaces; ++faceI)
    {
        const vector n = F.normals()[faceI];
        const label nGauss = F.nGaussPerFace();
        const List<scalar>& w = F.faces()[faceI]->weights();
        const List<scalar>& J2D = F.faces()[faceI]->J2D();

        for (label g = 0; g < nGauss; ++g)
        {
            // F(g) is vector ⇒ dot with n to get scalar flux
            surfInt += J2D[g] * phi.minusValueOnFace(faceI, g) * (F.fluxOnFace(faceI, g) & n) * w[g];
        }
    }

    return surfInt;
}

/**
 * \brief Surface integral for tensor flux
 * 
 * \param phi Scalar field at face Gauss points
 * \param F   Vector field at face Gauss points
 */
inline vector surfaceIntegral(const faceGaussField<scalar>& phi, const faceGaussField<tensor>& F)
{
    vector surfInt = Zero;

    const label nFaces = F.nFaces();

    for (label faceI = 0; faceI < nFaces; ++faceI)
    {
        const label nGauss = F.nGaussPerFace();
        const List<scalar>& w = F.faces()[faceI]->weights();
        const List<scalar>& J2D = F.faces()[faceI]->J2D();

        for (label g = 0; g < nGauss; ++g)
        {
            // fluxOnFace already returns T·n (vector)
            surfInt += J2D[g] * phi.minusValueOnFace(faceI, g) * F.fluxOnFace(faceI, g) * w[g];
        }
    }

    return surfInt;
}

/**
 * \brief Surface integral for gradient-projection terms using face mid-state.
 *
 * Computes \f$\int_{\partial\Omega}\phi\,C^\ast n\,dS\f$.
 * If C is scalar the result is vector; if C is vector the result is tensor.
 *
 * \param phi Scalar basis function at face Gauss points
 * \param C   Projected field at face Gauss points
 */
template<class Type>
inline typename outerProduct<Type, vector>::type gradProjectionSurfaceIntegral
(
    const faceGaussField<scalar>& phi,
    const faceGaussField<Type>& C
)
{
    typename outerProduct<Type, vector>::type surfInt =
        pTraits<typename outerProduct<Type, vector>::type>::zero;

    const label nFaces = C.nFaces();

    for (label faceI = 0; faceI < nFaces; ++faceI)
    {
        const vector n = C.normals()[faceI];
        const label nGauss = C.nGaussPerFace();
        const List<scalar>& w = C.faces()[faceI]->weights();
        const List<scalar>& J2D = C.faces()[faceI]->J2D();

        for (label g = 0; g < nGauss; ++g)
        {
            surfInt +=
                J2D[g]
              * phi.minusValueOnFace(faceI, g)
              * (C.midValueOnFace(faceI, g) * n)
              * w[g];
        }
    }

    return surfInt;
}

} // End namespace Foam

#endif

// ************************************************************************* //
