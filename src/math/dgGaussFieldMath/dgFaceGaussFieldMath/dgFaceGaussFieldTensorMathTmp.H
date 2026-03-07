/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Tensor operations for faceGaussField (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldTensorMathTmp_H
#define Foam_dgFaceGaussFieldTensorMathTmp_H

#include "faceGaussField.H"
#include "Tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "scalar.H"
#include "dgFaceGaussFieldTmpOpCore.H"

namespace Foam
{

struct DgFaceTensorEvalDev
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(dev(a)) { return dev(a); }
};

struct DgFaceTensorEvalSymm
{
    inline symmTensor operator()(const tensor& a) const { return symm(a); }
};

struct DgFaceTensorEvalSkew
{
    inline tensor operator()(const tensor& a) const { return skew(a); }
};

struct DgFaceTensorEvalTr
{
    template<class A>
    inline scalar operator()(const A& a) const { return tr(a); }
};

struct DgFaceTensorEvalDet
{
    inline scalar operator()(const tensor& a) const { return det(a); }
};

struct DgFaceTensorEvalT
{
    inline tensor operator()(const tensor& a) const { return T(a); }
};

struct DgFaceTensorEvalInv
{
    template<class A>
    inline tensor operator()(const A& a) const { return inv(a); }
};

struct DgFaceTensorMutDevTensor
{
    inline void operator()(tensor& ap, tensor& am) const
    {
        ap = dev(ap);
        am = dev(am);
    }
};

struct DgFaceTensorMutSkewTensor
{
    inline void operator()(tensor& ap, tensor& am) const
    {
        ap = skew(ap);
        am = skew(am);
    }
};

struct DgFaceTensorMutTTensor
{
    inline void operator()(tensor& ap, tensor& am) const
    {
        ap = T(ap);
        am = T(am);
    }
};

struct DgFaceTensorMutDevSymm
{
    inline void operator()(symmTensor& ap, symmTensor& am) const
    {
        ap = dev(ap);
        am = dev(am);
    }
};

struct DgFaceTensorEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

struct DgFaceTensorEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

struct DgFaceTensorEvalAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b)
    {
        return a & b;
    }
};

struct DgFaceTensorEvalAndRev
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(b & a)
    {
        return b & a;
    }
};

struct DgFaceTensorEvalAndAnd
{
    template<class A, class B>
    inline scalar operator()(const A& a, const B& b) const
    {
        return a && b;
    }
};

struct DgFaceTensorEvalAndAndRev
{
    template<class A, class B>
    inline scalar operator()(const A& a, const B& b) const
    {
        return b && a;
    }
};

struct DgFaceTensorMutAddA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap = ap + bp;
        am = am + bm;
    }
};

struct DgFaceTensorMutAddB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = ap + bp;
        bm = am + bm;
    }
};

struct DgFaceTensorMutSubA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap = ap - bp;
        am = am - bm;
    }
};

struct DgFaceTensorMutSubB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = ap - bp;
        bm = am - bm;
    }
};

struct DgFaceTensorMutAndA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap = ap & bp;
        am = am & bm;
    }
};

struct DgFaceTensorMutAndB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = ap & bp;
        bm = am & bm;
    }
};

struct DgFaceTensorMutAndBRev
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = bp & ap;
        bm = bm & am;
    }
};

// tensor utilities
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(dev,  faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalDev,  DgFaceTensorMutDevTensor)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(symm, faceGaussField<symmTensor>, faceGaussField<tensor>,     DgFaceTensorEvalSymm, Foam::dgFaceTmp::NoMutateUnary)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(skew, faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalSkew, DgFaceTensorMutSkewTensor)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(tr,   faceGaussField<scalar>,     faceGaussField<tensor>,     DgFaceTensorEvalTr,   Foam::dgFaceTmp::NoMutateUnary)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(det,  faceGaussField<scalar>,     faceGaussField<tensor>,     DgFaceTensorEvalDet,  Foam::dgFaceTmp::NoMutateUnary)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(T,    faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalT,    DgFaceTensorMutTTensor)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(inv,  faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalInv,  Foam::dgFaceTmp::NoMutateUnary)

// symmTensor utilities (no T(symmTensor) in original face file)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(tr,   faceGaussField<scalar>,     faceGaussField<symmTensor>, DgFaceTensorEvalTr,   Foam::dgFaceTmp::NoMutateUnary)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(dev,  faceGaussField<symmTensor>, faceGaussField<symmTensor>, DgFaceTensorEvalDev,  DgFaceTensorMutDevSymm)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION(inv,  faceGaussField<tensor>,     faceGaussField<symmTensor>, DgFaceTensorEvalInv,  Foam::dgFaceTmp::NoMutateUnary)

// tensor - tensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(+,  faceGaussField<tensor>, faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalAdd,    DgFaceTensorMutAddA, DgFaceTensorMutAddB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(-,  faceGaussField<tensor>, faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalSub,    DgFaceTensorMutSubA, DgFaceTensorMutSubB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&,  faceGaussField<tensor>, faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalAnd,    DgFaceTensorMutAndA, DgFaceTensorMutAndB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&&, faceGaussField<scalar>, faceGaussField<tensor>,     faceGaussField<tensor>,     DgFaceTensorEvalAndAnd, Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)

// tensor - symmTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(+,  faceGaussField<tensor>, faceGaussField<tensor>,     faceGaussField<symmTensor>, DgFaceTensorEvalAdd,    DgFaceTensorMutAddA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(-,  faceGaussField<tensor>, faceGaussField<tensor>,     faceGaussField<symmTensor>, DgFaceTensorEvalSub,    DgFaceTensorMutSubA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&,  faceGaussField<tensor>, faceGaussField<tensor>,     faceGaussField<symmTensor>, DgFaceTensorEvalAnd,    DgFaceTensorMutAndA, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&&, faceGaussField<scalar>, faceGaussField<tensor>,     faceGaussField<symmTensor>, DgFaceTensorEvalAndAnd, Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)

// symmTensor - tensor (preserve original semantics for &, && as B & A, B && A)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(+,  faceGaussField<tensor>, faceGaussField<symmTensor>, faceGaussField<tensor>,     DgFaceTensorEvalAdd,    Foam::dgFaceTmp::NoMutateLeft, DgFaceTensorMutAddB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(-,  faceGaussField<tensor>, faceGaussField<symmTensor>, faceGaussField<tensor>,     DgFaceTensorEvalSub,    Foam::dgFaceTmp::NoMutateLeft, DgFaceTensorMutSubB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&,  faceGaussField<tensor>, faceGaussField<symmTensor>, faceGaussField<tensor>,     DgFaceTensorEvalAndRev, Foam::dgFaceTmp::NoMutateLeft, DgFaceTensorMutAndBRev)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&&, faceGaussField<scalar>, faceGaussField<symmTensor>, faceGaussField<tensor>,     DgFaceTensorEvalAndAndRev, Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)

// symmTensor - symmTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(+,  faceGaussField<symmTensor>, faceGaussField<symmTensor>, faceGaussField<symmTensor>, DgFaceTensorEvalAdd,    DgFaceTensorMutAddA, DgFaceTensorMutAddB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(-,  faceGaussField<symmTensor>, faceGaussField<symmTensor>, faceGaussField<symmTensor>, DgFaceTensorEvalSub,    DgFaceTensorMutSubA, DgFaceTensorMutSubB)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&,  faceGaussField<tensor>,     faceGaussField<symmTensor>, faceGaussField<symmTensor>, DgFaceTensorEvalAnd,    Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)
DG_FACE_TMP_DEFINE_BINARY_OPERATOR(&&, faceGaussField<scalar>,     faceGaussField<symmTensor>, faceGaussField<symmTensor>, DgFaceTensorEvalAndAnd, Foam::dgFaceTmp::NoMutateLeft, Foam::dgFaceTmp::NoMutateRight)

} // End namespace Foam

#endif

// ************************************************************************* //

