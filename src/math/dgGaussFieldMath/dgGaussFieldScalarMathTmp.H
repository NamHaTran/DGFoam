/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

template<class OutType, class RefType>
inline tmp<GaussField<OutType>> dgGaussScalarNewLike(const GaussField<RefType>& A)
{
    return tmp<GaussField<OutType>>(new GaussField<OutType>(A.cellID(), A.dgMesh()));
}

struct DgGaussEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

struct DgGaussEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

struct DgGaussEvalMul
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b)
    {
        return a * b;
    }
};

struct DgGaussEvalDiv
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a / b)
    {
        return a / b;
    }
};

struct DgGaussMutAddLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a + b; }
};

struct DgGaussMutAddRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a + b; }
};

struct DgGaussMutSubLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a - b; }
};

struct DgGaussMutSubRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

struct DgGaussMutMulLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a * b; }
};

struct DgGaussMutMulRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a * b; }
};

struct DgGaussMutDivLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a / b; }
};

struct DgGaussMutDivRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a / b; }
};

struct DgGaussNoMutateLeft
{
    template<class A, class B>
    inline void operator()(A&, const B&) const {}
};

struct DgGaussNoMutateRight
{
    template<class A, class B>
    inline void operator()(const A&, B&) const {}
};

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

struct DgGaussPredLess
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a < b; }
};

struct DgGaussPredLessEq
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a <= b; }
};

struct DgGaussPredGreater
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a > b; }
};

struct DgGaussPredGreaterEq
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a >= b; }
};

struct DgGaussPredEq
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a == b; }
};

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

// Scalar - Scalar
DG_GAUSS_TMP_DEFINE_BINARY(+, scalar, scalar, scalar, DgGaussEvalAdd, DgGaussMutAddLeft, DgGaussMutAddRight)
DG_GAUSS_TMP_DEFINE_BINARY(-, scalar, scalar, scalar, DgGaussEvalSub, DgGaussMutSubLeft, DgGaussMutSubRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, scalar, scalar, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, scalar, scalar, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussMutDivRight)

// Scalar - Vector
DG_GAUSS_TMP_DEFINE_BINARY(*, vector, scalar, vector, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, vector, vector, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, vector, vector, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

// Scalar - Tensor
DG_GAUSS_TMP_DEFINE_BINARY(*, tensor, scalar, tensor, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, tensor, tensor, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, tensor, tensor, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

// Scalar - SymmTensor
DG_GAUSS_TMP_DEFINE_BINARY(*, symmTensor, scalar, symmTensor, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, symmTensor, symmTensor, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, symmTensor, symmTensor, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

// Scalar - SphericalTensor
DG_GAUSS_TMP_DEFINE_BINARY(*, sphericalTensor, scalar, sphericalTensor, DgGaussEvalMul, DgGaussNoMutateLeft, DgGaussMutMulRight)
DG_GAUSS_TMP_DEFINE_BINARY(*, sphericalTensor, sphericalTensor, scalar, DgGaussEvalMul, DgGaussMutMulLeft, DgGaussNoMutateRight)
DG_GAUSS_TMP_DEFINE_BINARY(/, sphericalTensor, sphericalTensor, scalar, DgGaussEvalDiv, DgGaussMutDivLeft, DgGaussNoMutateRight)

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
