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
 * \brief Volume integral for the sum of vector flux fields.
 *
 * Computes \f$\int_K \nabla\phi \cdot \sum_m F_m\,dV\f$ without building a
 * temporary summed field.
 *
 * \param gradPhi Gradient of basis function at cell Gauss points
 * \param fields  Pointer list of vector flux fields
 */
template<class FieldPtrList>
inline scalar volumeIntegral
(
    const cellGaussField<vector>& gradPhi,
    const FieldPtrList& fields,
    const label nFields
)
{
    scalar volInt = Zero;

    if (nFields == 0)
    {
        return volInt;
    }

    const label nGauss = gradPhi.size();
    const List<scalar>& w = gradPhi.cell()->weights();
    const List<scalar>& J3D = gradPhi.cell()->J3D();

    for (label g = 0; g < nGauss; ++g)
    {
        scalar fluxAtGp = Zero;

        for (label fieldI = 0; fieldI < nFields; ++fieldI)
        {
            fluxAtGp += gradPhi[g] & (*fields[fieldI])[g];
        }

        volInt += J3D[g] * fluxAtGp * w[g];
    }

    return volInt;
}

/**
 * \brief Volume integral for the sum of tensor flux fields.
 *
 * Computes \f$\int_K \nabla\phi \cdot \sum_m F_m\,dV\f$ without building a
 * temporary summed field.
 *
 * \param gradPhi Gradient of basis function at cell Gauss points
 * \param fields  Pointer list of tensor flux fields
 */
template<class FieldPtrList>
inline vector volumeIntegral
(
    const cellGaussField<vector>& gradPhi,
    const FieldPtrList& fields,
    const label nFields,
    const tensor*
)
{
    vector volInt = Zero;

    if (nFields == 0)
    {
        return volInt;
    }

    const label nGauss = gradPhi.size();
    const List<scalar>& w = gradPhi.cell()->weights();
    const List<scalar>& J3D = gradPhi.cell()->J3D();

    for (label g = 0; g < nGauss; ++g)
    {
        vector fluxAtGp = Zero;

        for (label fieldI = 0; fieldI < nFields; ++fieldI)
        {
            fluxAtGp += gradPhi[g] & (*fields[fieldI])[g];
        }

        volInt += J3D[g] * fluxAtGp * w[g];
    }

    return volInt;
}

/**
 * \brief Volume integrals for several gradient-projection state fields.
 *
 * Computes, for each active field \f$C_m\f$,
 *
 * \f[
 *     \int_K C_m \nabla\phi\,dV
 * \f]
 *
 * in a shared quadrature loop.  No temporary summed field is built since each
 * projection term produces a separate residual vector.
 *
 * \param gradPhi Gradient of basis function at cell Gauss points
 * \param fields  Pointer list of projected state fields
 * \param nFields Number of active entries in \p fields and \p results
 * \param results Output integrals, one per active field
 */
template<class Type, class FieldPtrList, class ResultList>
inline void gradProjectionVolumeIntegrals
(
    const cellGaussField<vector>& gradPhi,
    const FieldPtrList& fields,
    const label nFields,
    ResultList& results
)
{
    using result_type = typename outerProduct<Type, vector>::type;

    for (label fieldI = 0; fieldI < nFields; ++fieldI)
    {
        results[fieldI] = pTraits<result_type>::zero;
    }

    if (nFields == 0)
    {
        return;
    }

    const label nGauss = gradPhi.size();
    const List<scalar>& w = gradPhi.cell()->weights();
    const List<scalar>& J3D = gradPhi.cell()->J3D();

    for (label g = 0; g < nGauss; ++g)
    {
        const scalar Jw = J3D[g]*w[g];

        for (label fieldI = 0; fieldI < nFields; ++fieldI)
        {
            results[fieldI] += Jw*((*fields[fieldI])[g]*gradPhi[g]);
        }
    }
}

/**
 * \brief Surface integral for the sum of vector numerical flux fields.
 *
 * Computes \f$\int_{\partial K} \phi \sum_m F^\ast_m\cdot n\,dS\f$ without
 * building a temporary summed field.
 *
 * \param phi    Scalar basis function at face Gauss points
 * \param fields Pointer list of vector face flux fields
 */
template<class FieldPtrList>
inline scalar surfaceIntegral
(
    const faceGaussField<scalar>& phi,
    const FieldPtrList& fields,
    const label nFields
)
{
    scalar surfInt = Zero;

    if (nFields == 0)
    {
        return surfInt;
    }

    const faceGaussField<vector>& first = *fields[0];

    for (label faceI = 0; faceI < first.nFaces(); ++faceI)
    {
        const vector n = first.normals()[faceI];
        const label nGauss = first.nGaussPerFace();
        const List<scalar>& w = first.faces()[faceI]->weights();
        const List<scalar>& J2D = first.faces()[faceI]->J2D();

        for (label g = 0; g < nGauss; ++g)
        {
            scalar fluxAtFace = Zero;

            for (label fieldI = 0; fieldI < nFields; ++fieldI)
            {
                fluxAtFace += (*fields[fieldI]).fluxOnFace(faceI, g) & n;
            }

            surfInt +=
                J2D[g] * phi.minusValueOnFace(faceI, g) * fluxAtFace * w[g];
        }
    }

    return surfInt;
}

/**
 * \brief Surface integral for the sum of tensor numerical flux fields.
 *
 * For tensor fluxes, \c fluxOnFace() already stores
 * \f$F^\ast\cdot n\f$.
 *
 * \param phi    Scalar basis function at face Gauss points
 * \param fields Pointer list of tensor face flux fields
 */
template<class FieldPtrList>
inline vector surfaceIntegral
(
    const faceGaussField<scalar>& phi,
    const FieldPtrList& fields,
    const label nFields,
    const tensor*
)
{
    vector surfInt = Zero;

    if (nFields == 0)
    {
        return surfInt;
    }

    const faceGaussField<tensor>& first = *fields[0];

    for (label faceI = 0; faceI < first.nFaces(); ++faceI)
    {
        const label nGauss = first.nGaussPerFace();
        const List<scalar>& w = first.faces()[faceI]->weights();
        const List<scalar>& J2D = first.faces()[faceI]->J2D();

        for (label g = 0; g < nGauss; ++g)
        {
            vector fluxAtFace = Zero;

            for (label fieldI = 0; fieldI < nFields; ++fieldI)
            {
                fluxAtFace += (*fields[fieldI]).fluxOnFace(faceI, g);
            }

            surfInt +=
                J2D[g] * phi.minusValueOnFace(faceI, g) * fluxAtFace * w[g];
        }
    }

    return surfInt;
}

/**
 * \brief Volume integral for the sum of source fields.
 *
 * Computes \f$\int_K \phi \sum_m S_m\,dV\f$ without building a temporary
 * summed source field.
 *
 * \param phi    Scalar basis function at cell Gauss points
 * \param fields Pointer list of source fields
 */
template<class Type, class FieldPtrList>
inline Type sourceIntegral
(
    const cellGaussField<scalar>& phi,
    const FieldPtrList& fields,
    const label nFields
)
{
    Type sourceInt = pTraits<Type>::zero;

    if (nFields == 0)
    {
        return sourceInt;
    }

    const cellGaussField<Type>& first = *fields[0];
    const List<scalar>& w = first.cell()->weights();
    const List<scalar>& J3D = first.cell()->J3D();

    for (label g = 0; g < first.size(); ++g)
    {
        Type sourceAtGp = pTraits<Type>::zero;

        for (label fieldI = 0; fieldI < nFields; ++fieldI)
        {
            sourceAtGp += (*fields[fieldI])[g];
        }

        sourceInt += (J3D[g] * w[g] * phi[g]) * sourceAtGp;
    }

    return sourceInt;
}

/**
 * \brief Central surface integrals for several gradient-projection fields.
 *
 * Computes, for each active field \f$C_m\f$,
 *
 * \f[
 *     \int_{\partial K}\phi C_m^\ast n\,dS
 * \f]
 *
 * in a shared quadrature loop.  Boundary faces use the supplied middle state;
 * internal and processor faces use
 * \f$C^\ast = 0.5(C^- + C^+)\f$ on read.
 *
 * \param phi     Scalar basis function at face Gauss points
 * \param fields  Pointer list of projected face fields
 * \param nFields Number of active entries in \p fields and \p results
 * \param results Output integrals, one per active field
 */
template<class Type, class FieldPtrList, class ResultList>
inline void gradProjectionCentralSurfaceIntegrals
(
    const faceGaussField<scalar>& phi,
    const FieldPtrList& fields,
    const label nFields,
    ResultList& results
)
{
    using result_type = typename outerProduct<Type, vector>::type;

    for (label fieldI = 0; fieldI < nFields; ++fieldI)
    {
        results[fieldI] = pTraits<result_type>::zero;
    }

    if (nFields == 0)
    {
        return;
    }

    const faceGaussField<Type>& first = *fields[0];
    const label nFaces = first.nFaces();

    for (label faceI = 0; faceI < nFaces; ++faceI)
    {
        const vector n = first.normals()[faceI];
        const label nGauss = first.nGaussPerFace();
        const List<scalar>& w = first.faces()[faceI]->weights();
        const List<scalar>& J2D = first.faces()[faceI]->J2D();

        const bool useBoundaryMid =
            first.isBoundary(faceI) && !first.isProcessorPatch(faceI);

        for (label g = 0; g < nGauss; ++g)
        {
            const scalar JwPhi =
                J2D[g]*w[g]*phi.minusValueOnFace(faceI, g);

            for (label fieldI = 0; fieldI < nFields; ++fieldI)
            {
                const faceGaussField<Type>& C = *fields[fieldI];

                const Type CStar =
                    useBoundaryMid
                  ? C.midValueOnFace(faceI, g)
                  : scalar(0.5)
                   *(
                        C.minusValueOnFace(faceI, g)
                      + C.plusValueOnFace(faceI, g)
                    );

                results[fieldI] += JwPhi*(CStar*n);
            }
        }
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
