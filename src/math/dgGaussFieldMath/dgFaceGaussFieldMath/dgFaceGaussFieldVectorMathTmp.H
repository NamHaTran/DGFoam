/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

struct DgFaceVecEvalAdd
{
    inline vector operator()(const vector& a, const vector& b) const { return a + b; }
};

struct DgFaceVecEvalSub
{
    inline vector operator()(const vector& a, const vector& b) const { return a - b; }
};

struct DgFaceVecEvalCross
{
    inline vector operator()(const vector& a, const vector& b) const { return a ^ b; }
};

struct DgFaceVecEvalDot
{
    inline scalar operator()(const vector& a, const vector& b) const { return a & b; }
};

struct DgFaceVecEvalOuter
{
    inline tensor operator()(const vector& a, const vector& b) const { return a * b; }
};

struct DgFaceVecEvalAndVT
{
    template<class T>
    inline vector operator()(const vector& a, const T& b) const { return a & b; }
};

struct DgFaceVecEvalAndTV
{
    template<class T>
    inline vector operator()(const T& a, const vector& b) const { return a & b; }
};

struct DgFaceVecEvalMag
{
    inline scalar operator()(const vector& a) const { return mag(a); }
};

struct DgFaceVecEvalMagSqr
{
    inline scalar operator()(const vector& a) const { return magSqr(a); }
};

struct DgFaceVecMutAddA
{
    inline void operator()(vector& ap, vector& am, const vector& bp, const vector& bm) const
    {
        ap += bp;
        am += bm;
    }
};

struct DgFaceVecMutAddB
{
    inline void operator()(const vector& ap, const vector& am, vector& bp, vector& bm) const
    {
        bp += ap;
        bm += am;
    }
};

struct DgFaceVecMutSubA
{
    inline void operator()(vector& ap, vector& am, const vector& bp, const vector& bm) const
    {
        ap -= bp;
        am -= bm;
    }
};

struct DgFaceVecMutSubB
{
    inline void operator()(const vector& ap, const vector& am, vector& bp, vector& bm) const
    {
        bp = ap - bp;
        bm = am - bm;
    }
};

struct DgFaceVecMutCrossA
{
    inline void operator()(vector& ap, vector& am, const vector& bp, const vector& bm) const
    {
        ap = ap ^ bp;
        am = am ^ bm;
    }
};

struct DgFaceVecMutCrossB
{
    inline void operator()(const vector& ap, const vector& am, vector& bp, vector& bm) const
    {
        bp = ap ^ bp;
        bm = am ^ bm;
    }
};

struct DgFaceVecMutAndA
{
    template<class T>
    inline void operator()(vector& ap, vector& am, const T& bp, const T& bm) const
    {
        ap = ap & bp;
        am = am & bm;
    }
};

struct DgFaceVecMutAndB
{
    template<class T>
    inline void operator()(const T& ap, const T& am, vector& bp, vector& bm) const
    {
        bp = ap & bp;
        bm = am & bm;
    }
};

// Vector - Vector
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(+, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalAdd,   DgFaceVecMutAddA,   DgFaceVecMutAddB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(-, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalSub,   DgFaceVecMutSubA,   DgFaceVecMutSubB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(^, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalCross, DgFaceVecMutCrossA, DgFaceVecMutCrossB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<scalar>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalDot,   Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(*, faceGaussField<tensor>, faceGaussField<vector>, faceGaussField<vector>, DgFaceVecEvalOuter, Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)

DG_FACE_TMP_DEFINE_UNARY_FUNCTION(mag,    faceGaussField<scalar>, faceGaussField<vector>, DgFaceVecEvalMag,    Foam::dgFaceTmp::NoMutateUnary)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(magSqr, faceGaussField<scalar>, faceGaussField<vector>, DgFaceVecEvalMagSqr, Foam::dgFaceTmp::NoMutateUnary)

// Vector - Tensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<tensor>, DgFaceVecEvalAndVT, DgFaceVecMutAndA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<tensor>, faceGaussField<vector>, DgFaceVecEvalAndTV, Foam::dgFaceTmp::NoMutateLeft, DgFaceVecMutAndB)

// Vector - symmTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<symmTensor>, DgFaceVecEvalAndVT, DgFaceVecMutAndA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<symmTensor>, faceGaussField<vector>, DgFaceVecEvalAndTV, Foam::dgFaceTmp::NoMutateLeft, DgFaceVecMutAndB)

// Vector - sphericalTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<vector>, faceGaussField<sphericalTensor>, DgFaceVecEvalAndVT, DgFaceVecMutAndA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&, faceGaussField<vector>, faceGaussField<sphericalTensor>, faceGaussField<vector>, DgFaceVecEvalAndTV, Foam::dgFaceTmp::NoMutateLeft, DgFaceVecMutAndB)

} // End namespace Foam

#endif

// ************************************************************************* //

