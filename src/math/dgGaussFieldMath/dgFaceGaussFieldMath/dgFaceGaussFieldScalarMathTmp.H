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

/** \brief Evaluate binary addition at one face Gauss point. */
struct DgFaceEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

/** \brief Evaluate binary subtraction at one face Gauss point. */
struct DgFaceEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

/** \brief Evaluate binary multiplication at one face Gauss point. */
struct DgFaceEvalMul
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b)
    {
        return a * b;
    }
};

/** \brief Evaluate binary division at one face Gauss point. */
struct DgFaceEvalDiv
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a / b)
    {
        return a / b;
    }
};

/** \brief In-place mutate left operand for addition. */
struct DgFaceMutAddA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap += bp;
        am += bm;
    }
};

/** \brief In-place mutate right operand for addition. */
struct DgFaceMutAddB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp += ap;
        bm += am;
    }
};

/** \brief In-place mutate left operand for subtraction. */
struct DgFaceMutSubA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap -= bp;
        am -= bm;
    }
};

/** \brief In-place mutate right operand for subtraction. */
struct DgFaceMutSubB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = ap - bp;
        bm = am - bm;
    }
};

/** \brief In-place mutate left operand for multiplication. */
struct DgFaceMutMulA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap *= bp;
        am *= bm;
    }
};

/** \brief In-place mutate right operand for multiplication. */
struct DgFaceMutMulB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp *= ap;
        bm *= am;
    }
};

/** \brief In-place mutate left operand for division. */
struct DgFaceMutDivA
{
    template<class A, class B>
    inline void operator()(A& ap, A& am, const B& bp, const B& bm) const
    {
        ap /= bp;
        am /= bm;
    }
};

/** \brief In-place mutate right operand for division. */
struct DgFaceMutDivB
{
    template<class A, class B>
    inline void operator()(const A& ap, const A& am, B& bp, B& bm) const
    {
        bp = ap / bp;
        bm = am / bm;
    }
};

/** \brief Predicate functor for `<` on scalar face fields. */
struct DgFacePredLess
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

/** \brief Predicate functor for `<=` on scalar face fields. */
struct DgFacePredLessEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

/** \brief Predicate functor for `>` on scalar face fields. */
struct DgFacePredGreater
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

/** \brief Predicate functor for `>=` on scalar face fields. */
struct DgFacePredGreaterEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

/** \brief Predicate functor for `==` on scalar face fields. */
struct DgFacePredEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

/** \brief Predicate functor for `!=` on scalar face fields. */
struct DgFacePredNe
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

/** \brief Scalar-scalar operators: `+`, `-`, `*`, `/`. */
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
