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
    Operator overloads between GaussField<T> types (tmp-aware):
    - scalar - scalar
    - scalar - vector
    - scalar - tensor
    - scalar - symmTensor
    - scalar - sphericalTensor

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldScalarMathTmp_H
#define Foam_dgGaussFieldScalarMathTmp_H

#include "GaussField.H"
#include "dgFaceGaussFieldScalarMathTmp.H"
#include "dgCellGaussFieldScalarMathTmp.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include <type_traits>

namespace Foam
{

inline scalar maxValue(const tmp<GaussField<scalar>>& tA)
{
    const GaussField<scalar>& A = tA();
    scalar result = max(maxValue(A.cellField()), maxValue(A.faceField()));
    tA.clear();
    return result;
}

inline scalar minValue(const tmp<GaussField<scalar>>& tA)
{
    const GaussField<scalar>& A = tA();
    scalar result = min(minValue(A.cellField()), minValue(A.faceField()));
    tA.clear();
    return result;
}

/** \brief Allocate output GaussField with geometry copied from reference field. */
template<class OutType, class RefType>
inline tmp<GaussField<OutType>> dgGaussScalarNewLike(const GaussField<RefType>& A)
{
    return tmp<GaussField<OutType>>(new GaussField<OutType>(A.cellID(), A.dgMesh()));
}

/** \brief Evaluate binary addition at one Gauss point. */
struct DgGaussEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

/** \brief Evaluate binary subtraction at one Gauss point. */
struct DgGaussEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

/** \brief Evaluate binary multiplication at one Gauss point. */
struct DgGaussEvalMul
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b)
    {
        return a * b;
    }
};

/** \brief Evaluate binary division at one Gauss point. */
struct DgGaussEvalDiv
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a / b)
    {
        return a / b;
    }
};

/** \brief In-place mutate left operand for addition. */
struct DgGaussMutAddLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a + b; }
};

/** \brief In-place mutate right operand for addition. */
struct DgGaussMutAddRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a + b; }
};

/** \brief In-place mutate left operand for subtraction. */
struct DgGaussMutSubLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a - b; }
};

/** \brief In-place mutate right operand for subtraction. */
struct DgGaussMutSubRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

/** \brief In-place mutate left operand for multiplication. */
struct DgGaussMutMulLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a * b; }
};

/** \brief In-place mutate right operand for multiplication. */
struct DgGaussMutMulRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a * b; }
};

/** \brief In-place mutate left operand for division. */
struct DgGaussMutDivLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a / b; }
};

/** \brief In-place mutate right operand for division. */
struct DgGaussMutDivRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a / b; }
};

/** \brief No-op mutator for left operand when tmp reuse is disabled. */
struct DgGaussNoMutateLeft
{
    template<class A, class B>
    inline void operator()(A&, const B&) const {}
};

/** \brief No-op mutator for right operand when tmp reuse is disabled. */
struct DgGaussNoMutateRight
{
    template<class A, class B>
    inline void operator()(const A&, B&) const {}
};

/** \brief Core binary dispatcher for non-tmp GaussField inputs. */
template<class OutType, class AType, class BType, class Eval>
inline tmp<GaussField<OutType>> dgGaussBinary
(
    const GaussField<AType>& A,
    const GaussField<BType>& B,
    const Eval& eval
)
{
    auto tRes = dgGaussScalarNewLike<OutType>(A);
    GaussField<OutType>& R = tRes.ref();
    R.cellField() = eval(A.cellField(), B.cellField());
    R.faceField() = eval(A.faceField(), B.faceField());
    return tRes;
}

template<class OutType, class AType, class BType, class Eval, class MutA>
inline typename std::enable_if<std::is_same<OutType, AType>::value, tmp<GaussField<OutType>>>::type
dgGaussBinary
(
    const tmp<GaussField<AType>>& tA,
    const GaussField<BType>& B,
    const Eval& eval,
    const MutA& mutA
)
{
    if (tA.isTmp())
    {
        GaussField<AType>& A = const_cast<tmp<GaussField<AType>>&>(tA).ref();
        mutA(A.cellField(), B.cellField());
        mutA(A.faceField(), B.faceField());
        return tA;
    }

    return dgGaussBinary<OutType>(tA(), B, eval);
}

template<class OutType, class AType, class BType, class Eval, class MutA>
inline typename std::enable_if<!std::is_same<OutType, AType>::value, tmp<GaussField<OutType>>>::type
dgGaussBinary
(
    const tmp<GaussField<AType>>& tA,
    const GaussField<BType>& B,
    const Eval& eval,
    const MutA&
)
{
    return dgGaussBinary<OutType>(tA(), B, eval);
}

template<class OutType, class AType, class BType, class Eval, class MutB>
inline typename std::enable_if<std::is_same<OutType, BType>::value, tmp<GaussField<OutType>>>::type
dgGaussBinary
(
    const GaussField<AType>& A,
    const tmp<GaussField<BType>>& tB,
    const Eval& eval,
    const MutB& mutB
)
{
    if (tB.isTmp())
    {
        GaussField<BType>& B = const_cast<tmp<GaussField<BType>>&>(tB).ref();
        mutB(A.cellField(), B.cellField());
        mutB(A.faceField(), B.faceField());
        return tB;
    }

    return dgGaussBinary<OutType>(A, tB(), eval);
}

template<class OutType, class AType, class BType, class Eval, class MutB>
inline typename std::enable_if<!std::is_same<OutType, BType>::value, tmp<GaussField<OutType>>>::type
dgGaussBinary
(
    const GaussField<AType>& A,
    const tmp<GaussField<BType>>& tB,
    const Eval& eval,
    const MutB&
)
{
    return dgGaussBinary<OutType>(A, tB(), eval);
}

template<class OutType, class AType, class BType, class Eval, class MutA, class MutB>
inline tmp<GaussField<OutType>> dgGaussBinary
(
    const tmp<GaussField<AType>>& tA,
    const tmp<GaussField<BType>>& tB,
    const Eval& eval,
    const MutA& mutA,
    const MutB& mutB
)
{
    if (tA.isTmp())
    {
        return dgGaussBinary<OutType>(tA, tB(), eval, mutA);
    }
    return dgGaussBinary<OutType>(tA(), tB, eval, mutB);
}

template<class AType, class BType, class Pred>
inline bool dgGaussCompare
(
    const GaussField<AType>& A,
    const GaussField<BType>& B,
    const Pred& pred
)
{
    return pred(A.cellField(), B.cellField()) && pred(A.faceField(), B.faceField());
}

template<class AType, class BType, class Pred>
inline bool dgGaussCompare
(
    const tmp<GaussField<AType>>& tA,
    const GaussField<BType>& B,
    const Pred& pred
)
{
    return dgGaussCompare(tA(), B, pred);
}

template<class AType, class BType, class Pred>
inline bool dgGaussCompare
(
    const GaussField<AType>& A,
    const tmp<GaussField<BType>>& tB,
    const Pred& pred
)
{
    return dgGaussCompare(A, tB(), pred);
}

template<class AType, class BType, class Pred>
inline bool dgGaussCompare
(
    const tmp<GaussField<AType>>& tA,
    const tmp<GaussField<BType>>& tB,
    const Pred& pred
)
{
    return dgGaussCompare(tA(), tB(), pred);
}

/** \brief Predicate functor for `<` comparison. */
struct DgGaussPredLess
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a < b; }
};

/** \brief Predicate functor for `<=` comparison. */
struct DgGaussPredLessEq
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a <= b; }
};

/** \brief Predicate functor for `>` comparison. */
struct DgGaussPredGreater
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a > b; }
};

/** \brief Predicate functor for `>=` comparison. */
struct DgGaussPredGreaterEq
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a >= b; }
};

/** \brief Predicate functor for `==` comparison. */
struct DgGaussPredEq
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a == b; }
};

/** \brief Predicate functor for `!=` comparison. */
struct DgGaussPredNe
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a != b; }
};

#define DG_GAUSS_TMP_DEFINE_BINARY(OP, OUTF, AF, BF, EVAL, MUTA, MUTB)                  \
inline tmp<GaussField<OUTF>> operator OP (const GaussField<AF>& A, const GaussField<BF>& B) \
{                                                                                          \
    return dgGaussBinary<OUTF>(A, B, EVAL());                                             \
}                                                                                          \
inline tmp<GaussField<OUTF>> operator OP (const tmp<GaussField<AF>>& tA, const GaussField<BF>& B) \
{                                                                                          \
    return dgGaussBinary<OUTF>(tA, B, EVAL(), MUTA());                                    \
}                                                                                          \
inline tmp<GaussField<OUTF>> operator OP (const GaussField<AF>& A, const tmp<GaussField<BF>>& tB) \
{                                                                                          \
    return dgGaussBinary<OUTF>(A, tB, EVAL(), MUTB());                                    \
}                                                                                          \
inline tmp<GaussField<OUTF>> operator OP (const tmp<GaussField<AF>>& tA, const tmp<GaussField<BF>>& tB) \
{                                                                                          \
    return dgGaussBinary<OUTF>(tA, tB, EVAL(), MUTA(), MUTB());                           \
}

#define DG_GAUSS_TMP_DEFINE_COMPARE(OP, AF, BF, PRED)                                     \
inline bool operator OP (const GaussField<AF>& A, const GaussField<BF>& B)                \
{                                                                                           \
    return dgGaussCompare(A, B, PRED());                                                   \
}                                                                                           \
inline bool operator OP (const tmp<GaussField<AF>>& tA, const GaussField<BF>& B)          \
{                                                                                           \
    return dgGaussCompare(tA, B, PRED());                                                  \
}                                                                                           \
inline bool operator OP (const GaussField<AF>& A, const tmp<GaussField<BF>>& tB)          \
{                                                                                           \
    return dgGaussCompare(A, tB, PRED());                                                  \
}                                                                                           \
inline bool operator OP (const tmp<GaussField<AF>>& tA, const tmp<GaussField<BF>>& tB)    \
{                                                                                           \
    return dgGaussCompare(tA, tB, PRED());                                                 \
}

/** \brief Scalar-scalar operators: `+`, `-`, `*`, `/`. */
// Scalar - Scalar
DG_GAUSS_TMP_DEFINE_BINARY(+, scalar, scalar, scalar, DgGaussEvalAdd, DgGaussMutAddLeft, DgGaussMutAddRight)
DG_GAUSS_TMP_DEFINE_BINARY(-, scalar, scalar, scalar, DgGaussEvalSub, DgGaussMutSubLeft, DgGaussMutSubRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, scalar, scalar, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, scalar, scalar, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussMutDivRight)

/** \brief Scalar-vector operators: `*` and vector `/` scalar. */
// Scalar - Vector
DG_GAUSS_TMP_DEFINE_BINARY(*, vector, scalar, vector, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, vector, vector, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, vector, vector, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

/** \brief Scalar-tensor operators: `*` and tensor `/` scalar. */
// Scalar - Tensor
DG_GAUSS_TMP_DEFINE_BINARY(*, tensor, scalar, tensor, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, tensor, tensor, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, tensor, tensor, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

/** \brief Scalar-symmTensor operators: `*` and symmTensor `/` scalar. */
// Scalar - SymmTensor
DG_GAUSS_TMP_DEFINE_BINARY(*, symmTensor, scalar, symmTensor, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, symmTensor, symmTensor, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, symmTensor, symmTensor, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

/** \brief Scalar-sphericalTensor operators: `*` and sphericalTensor `/` scalar. */
// Scalar - SphericalTensor
DG_GAUSS_TMP_DEFINE_BINARY(*, sphericalTensor, scalar, sphericalTensor, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, sphericalTensor, sphericalTensor, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, sphericalTensor, sphericalTensor, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

/** \brief Comparison operators on scalar GaussField pairs. */
// Comparisons on scalar fields
DG_GAUSS_TMP_DEFINE_COMPARE(<,  scalar, scalar, DgGaussPredLess)
DG_GAUSS_TMP_DEFINE_COMPARE(<=, scalar, scalar, DgGaussPredLessEq)
DG_GAUSS_TMP_DEFINE_COMPARE(>,  scalar, scalar, DgGaussPredGreater)
DG_GAUSS_TMP_DEFINE_COMPARE(>=, scalar, scalar, DgGaussPredGreaterEq)
DG_GAUSS_TMP_DEFINE_COMPARE(==, scalar, scalar, DgGaussPredEq)
DG_GAUSS_TMP_DEFINE_COMPARE(!=, scalar, scalar, DgGaussPredNe)

#undef DG_GAUSS_TMP_DEFINE_BINARY
#undef DG_GAUSS_TMP_DEFINE_COMPARE

} // End namespace Foam

#endif

// ************************************************************************* //
