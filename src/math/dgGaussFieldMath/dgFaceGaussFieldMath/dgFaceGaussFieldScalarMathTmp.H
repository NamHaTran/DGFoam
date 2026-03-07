/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Operator overloads between faceGaussField<T> types (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldScalarMathTmp_H
#define Foam_dgFaceGaussFieldScalarMathTmp_H

#include "faceGaussField.H"
#include "dgFaceGaussFieldTmpOpCore.H"

namespace Foam
{

struct DgFaceEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

struct DgFaceEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

struct DgFaceEvalMul
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b)
    {
        return a * b;
    }
};

struct DgFaceEvalDiv
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a / b)
    {
        return a / b;
    }
};

struct DgFaceMutAddA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap += bp;
        am += bm;
    }
};

struct DgFaceMutAddB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp += ap;
        bm += am;
    }
};

struct DgFaceMutSubA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap -= bp;
        am -= bm;
    }
};

struct DgFaceMutSubB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = ap - bp;
        bm = am - bm;
    }
};

struct DgFaceMutMulA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap *= bp;
        am *= bm;
    }
};

struct DgFaceMutMulB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp *= ap;
        bm *= am;
    }
};

struct DgFaceMutDivA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap /= bp;
        am /= bm;
    }
};

struct DgFaceMutDivB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = ap / bp;
        bm = am / bm;
    }
};

struct DgFacePredLess
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

struct DgFacePredLessEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

struct DgFacePredGreater
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

struct DgFacePredGreaterEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

struct DgFacePredEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

struct DgFacePredNe
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

// Scalar - Scalar
DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    +,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFaceEvalAdd,
    DgFaceMutAddA,
    DgFaceMutAddB
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    -,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFaceEvalSub,
    DgFaceMutSubA,
    DgFaceMutSubB
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFaceEvalMul,
    DgFaceMutMulA,
    DgFaceMutMulB
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFaceEvalDiv,
    DgFaceMutDivA,
    DgFaceMutDivB
)

// Scalar - Vector
DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<vector>,
    faceGaussField<scalar>,
    faceGaussField<vector>,
    DgFaceEvalMul,
    Foam::dgFaceTmp::NoMutateLeft,
    DgFaceMutMulB
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<vector>,
    faceGaussField<vector>,
    faceGaussField<scalar>,
    DgFaceEvalMul,
    DgFaceMutMulA,
    Foam::dgFaceTmp::NoMutateRight
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    faceGaussField<vector>,
    faceGaussField<vector>,
    faceGaussField<scalar>,
    DgFaceEvalDiv,
    DgFaceMutDivA,
    Foam::dgFaceTmp::NoMutateRight
)

// Scalar - Tensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<tensor>,
    faceGaussField<scalar>,
    faceGaussField<tensor>,
    DgFaceEvalMul,
    Foam::dgFaceTmp::NoMutateLeft,
    DgFaceMutMulB
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<tensor>,
    faceGaussField<tensor>,
    faceGaussField<scalar>,
    DgFaceEvalMul,
    DgFaceMutMulA,
    Foam::dgFaceTmp::NoMutateRight
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    faceGaussField<tensor>,
    faceGaussField<tensor>,
    faceGaussField<scalar>,
    DgFaceEvalDiv,
    DgFaceMutDivA,
    Foam::dgFaceTmp::NoMutateRight
)

// Scalar - symmTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<symmTensor>,
    faceGaussField<scalar>,
    faceGaussField<symmTensor>,
    DgFaceEvalMul,
    Foam::dgFaceTmp::NoMutateLeft,
    DgFaceMutMulB
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<symmTensor>,
    faceGaussField<symmTensor>,
    faceGaussField<scalar>,
    DgFaceEvalMul,
    DgFaceMutMulA,
    Foam::dgFaceTmp::NoMutateRight
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    faceGaussField<symmTensor>,
    faceGaussField<symmTensor>,
    faceGaussField<scalar>,
    DgFaceEvalDiv,
    DgFaceMutDivA,
    Foam::dgFaceTmp::NoMutateRight
)

// Scalar - sphericalTensor
DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<sphericalTensor>,
    faceGaussField<scalar>,
    faceGaussField<sphericalTensor>,
    DgFaceEvalMul,
    Foam::dgFaceTmp::NoMutateLeft,
    DgFaceMutMulB
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    faceGaussField<sphericalTensor>,
    faceGaussField<sphericalTensor>,
    faceGaussField<scalar>,
    DgFaceEvalMul,
    DgFaceMutMulA,
    Foam::dgFaceTmp::NoMutateRight
)

DG_FACE_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    faceGaussField<sphericalTensor>,
    faceGaussField<sphericalTensor>,
    faceGaussField<scalar>,
    DgFaceEvalDiv,
    DgFaceMutDivA,
    Foam::dgFaceTmp::NoMutateRight
)

// Comparison operators between faceGaussField<scalar>
DG_FACE_TMP_DEFINE_COMPARE_OPERATOR
(
    <,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFacePredLess
)

DG_FACE_TMP_DEFINE_COMPARE_OPERATOR
(
    <=,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFacePredLessEq
)

DG_FACE_TMP_DEFINE_COMPARE_OPERATOR
(
    >,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFacePredGreater
)

DG_FACE_TMP_DEFINE_COMPARE_OPERATOR
(
    >=,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFacePredGreaterEq
)

DG_FACE_TMP_DEFINE_COMPARE_OPERATOR
(
    ==,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFacePredEq
)

DG_FACE_TMP_DEFINE_COMPARE_OPERATOR
(
    !=,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFacePredNe
)

} // End namespace Foam

#endif

// ************************************************************************* //

