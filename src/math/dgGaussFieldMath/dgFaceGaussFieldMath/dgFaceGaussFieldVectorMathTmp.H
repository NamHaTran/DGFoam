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

Description
    Vector math operations for faceGaussField (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldVectorMathTmp_H
#define Foam_dgFaceGaussFieldVectorMathTmp_H

#include "faceGaussField.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "VectorSpace.H"
#include "dgFaceGaussFieldTmpOpCore.H"

namespace Foam
{

/** \brief Evaluate vector addition at one face Gauss point. */
struct DgFaceVecEvalAdd
{
    inline vector operator()(const vector& a, const vector& b) const { return a + b; }
};

/** \brief Evaluate vector subtraction at one face Gauss point. */
struct DgFaceVecEvalSub
{
    inline vector operator()(const vector& a, const vector& b) const { return a - b; }
};

/** \brief Evaluate vector cross product at one face Gauss point. */
struct DgFaceVecEvalCross
{
    inline vector operator()(const vector& a, const vector& b) const { return a ^ b; }
};

/** \brief Evaluate vector dot product at one face Gauss point. */
struct DgFaceVecEvalDot
{
    inline scalar operator()(const vector& a, const vector& b) const { return a & b; }
};

/** \brief Evaluate vector outer product at one face Gauss point. */
struct DgFaceVecEvalOuter
{
    inline tensor operator()(const vector& a, const vector& b) const { return a * b; }
};

/** \brief Evaluate `vector & tensorLike`. */
struct DgFaceVecEvalAndVT
{
    template<class T>
    inline vector operator()(const vector& a, const T& b) const { return a & b; }
};

/** \brief Evaluate `tensorLike & vector`. */
struct DgFaceVecEvalAndTV
{
    template<class T>
    inline vector operator()(const T& a, const vector& b) const { return a & b; }
};

/** \brief Evaluate vector magnitude at one face Gauss point. */
struct DgFaceVecEvalMag
{
    inline scalar operator()(const vector& a) const { return mag(a); }
};

/** \brief Evaluate squared vector magnitude at one face Gauss point. */
struct DgFaceVecEvalMagSqr
{
    inline scalar operator()(const vector& a) const { return magSqr(a); }
};

/** \brief In-place mutate left operand for vector addition. */
struct DgFaceVecMutAddA
{
    inline void operator()(vector& ap, vector& am, const vector& bp, const vector& bm) const
    {
        ap += bp;
        am += bm;
    }
};

/** \brief In-place mutate right operand for vector addition. */
struct DgFaceVecMutAddB
{
    inline void operator()(const vector& ap, const vector& am, vector& bp, vector& bm) const
    {
        bp += ap;
        bm += am;
    }
};

/** \brief In-place mutate left operand for vector subtraction. */
struct DgFaceVecMutSubA
{
    inline void operator()(vector& ap, vector& am, const vector& bp, const vector& bm) const
    {
        ap -= bp;
        am -= bm;
    }
};

/** \brief In-place mutate right operand for vector subtraction. */
struct DgFaceVecMutSubB
{
    inline void operator()(const vector& ap, const vector& am, vector& bp, vector& bm) const
    {
        bp = ap - bp;
        bm = am - bm;
    }
};

/** \brief In-place mutate left operand for vector cross product. */
struct DgFaceVecMutCrossA
{
    inline void operator()(vector& ap, vector& am, const vector& bp, const vector& bm) const
    {
        ap = ap ^ bp;
        am = am ^ bm;
    }
};

/** \brief In-place mutate right operand for vector cross product. */
struct DgFaceVecMutCrossB
{
    inline void operator()(const vector& ap, const vector& am, vector& bp, vector& bm) const
    {
        bp = ap ^ bp;
        bm = am ^ bm;
    }
};

/** \brief In-place mutate left operand for `&` with tensor-like types. */
struct DgFaceVecMutAndA
{
    template<class T>
    inline void operator()(vector& ap, vector& am, const T& bp, const T& bm) const
    {
        ap = ap & bp;
        am = am & bm;
    }
};

/** \brief In-place mutate right operand for `&` with tensor-like types. */
struct DgFaceVecMutAndB
{
    template<class T>
    inline void operator()(const T& ap, const T& am, vector& bp, vector& bm) const
    {
        bp = ap & bp;
        bm = am & bm;
    }
};

/** \brief Binary vector-vector operators: `+`, `-`, `^`, `&`, `*`. */
// Vector - Vector
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(+, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalAdd,   DgFaceVecMutAddA,   DgFaceVecMutAddB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(-, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalSub,   DgFaceVecMutSubA,   DgFaceVecMutSubB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(^, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalCross, DgFaceVecMutCrossA, DgFaceVecMutCrossB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<scalar>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalDot,   Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(*, faceGaussField<tensor>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalOuter, Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)

/** \brief Unary vector operators: `mag`, `magSqr`. */
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(mag,    faceGaussField<scalar>, faceGaussField<vector>, DgFaceVecEvalMag,    Foam::dgFaceTmp::NoMutateUnary)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(magSqr, faceGaussField<scalar>, faceGaussField<vector>, DgFaceVecEvalMagSqr, Foam::dgFaceTmp::NoMutateUnary)

/** \brief Binary operators for `vector & tensor`. */
// Vector - Tensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<tensor>, DgFaceVecEvalAndVT, DgFaceVecMutAndA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<tensor>, faceGaussField<vector>, DgFaceVecEvalAndTV, Foam::dgFaceTmp::NoMutateLeft, DgFaceVecMutAndB)

/** \brief Binary operators for `vector & symmTensor`. */
// Vector - symmTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<symmTensor>, DgFaceVecEvalAndVT, DgFaceVecMutAndA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<symmTensor>, faceGaussField<vector>, DgFaceVecEvalAndTV, Foam::dgFaceTmp::NoMutateLeft, DgFaceVecMutAndB)

/** \brief Binary operators for `vector & sphericalTensor`. */
// Vector - sphericalTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<sphericalTensor>, DgFaceVecEvalAndVT, DgFaceVecMutAndA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<sphericalTensor>, faceGaussField<vector>, DgFaceVecEvalAndTV, Foam::dgFaceTmp::NoMutateLeft, DgFaceVecMutAndB)

} // End namespace Foam

#endif

// ************************************************************************* //
