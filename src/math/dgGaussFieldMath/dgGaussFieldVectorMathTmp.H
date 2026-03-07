/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Vector math operations for GaussField (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldVectorMathTmp_H
#define Foam_dgGaussFieldVectorMathTmp_H

#include "GaussField.H"
#include "dgFaceGaussFieldVectorMathTmp.H"
#include "dgCellGaussFieldVectorMathTmp.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include <type_traits>

namespace Foam
{

template<class OutType, class RefType>
inline tmp<GaussField<OutType>> dgGaussVectorNewLike(const GaussField<RefType>& A)
{
    return tmp<GaussField<OutType>>(new GaussField<OutType>(A.cellID(), A.dgMesh()));
}

template<class OutType, class InType, class Eval>
inline tmp<GaussField<OutType>> dgGaussVectorUnary
(
    const GaussField<InType>& A,
    const Eval& eval
)
{
    auto tRes = dgGaussVectorNewLike<OutType>(A);
    GaussField<OutType>& R = tRes.ref();
    R.cellField() = eval(A.cellField());
    R.faceField() = eval(A.faceField());
    return tRes;
}

template<class OutType, class InType, class Eval, class Mutate>
inline typename std::enable_if<std::is_same<OutType, InType>::value, tmp<GaussField<OutType>>>::type
dgGaussVectorUnary
(
    const tmp<GaussField<InType>>& tA,
    const Eval& eval,
    const Mutate& mutate
)
{
    if (tA.isTmp())
    {
        GaussField<InType>& A = const_cast<tmp<GaussField<InType>>&>(tA).ref();
        mutate(A.cellField());
        mutate(A.faceField());
        return tA;
    }
    return dgGaussVectorUnary<OutType>(tA(), eval);
}

template<class OutType, class InType, class Eval, class Mutate>
inline typename std::enable_if<!std::is_same<OutType, InType>::value, tmp<GaussField<OutType>>>::type
dgGaussVectorUnary
(
    const tmp<GaussField<InType>>& tA,
    const Eval& eval,
    const Mutate&
)
{
    return dgGaussVectorUnary<OutType>(tA(), eval);
}

template<class OutType, class AType, class BType, class Eval>
inline tmp<GaussField<OutType>> dgGaussVectorBinary
(
    const GaussField<AType>& A,
    const GaussField<BType>& B,
    const Eval& eval
)
{
    auto tRes = dgGaussVectorNewLike<OutType>(A);
    GaussField<OutType>& R = tRes.ref();
    R.cellField() = eval(A.cellField(), B.cellField());
    R.faceField() = eval(A.faceField(), B.faceField());
    return tRes;
}

template<class OutType, class AType, class BType, class Eval, class MutA>
inline typename std::enable_if<std::is_same<OutType, AType>::value, tmp<GaussField<OutType>>>::type
dgGaussVectorBinary
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
    return dgGaussVectorBinary<OutType>(tA(), B, eval);
}

template<class OutType, class AType, class BType, class Eval, class MutA>
inline typename std::enable_if<!std::is_same<OutType, AType>::value, tmp<GaussField<OutType>>>::type
dgGaussVectorBinary
(
    const tmp<GaussField<AType>>& tA,
    const GaussField<BType>& B,
    const Eval& eval,
    const MutA&
)
{
    return dgGaussVectorBinary<OutType>(tA(), B, eval);
}

template<class OutType, class AType, class BType, class Eval, class MutB>
inline typename std::enable_if<std::is_same<OutType, BType>::value, tmp<GaussField<OutType>>>::type
dgGaussVectorBinary
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
    return dgGaussVectorBinary<OutType>(A, tB(), eval);
}

template<class OutType, class AType, class BType, class Eval, class MutB>
inline typename std::enable_if<!std::is_same<OutType, BType>::value, tmp<GaussField<OutType>>>::type
dgGaussVectorBinary
(
    const GaussField<AType>& A,
    const tmp<GaussField<BType>>& tB,
    const Eval& eval,
    const MutB&
)
{
    return dgGaussVectorBinary<OutType>(A, tB(), eval);
}

template<class OutType, class AType, class BType, class Eval, class MutA, class MutB>
inline tmp<GaussField<OutType>> dgGaussVectorBinary
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
        return dgGaussVectorBinary<OutType>(tA, tB(), eval, mutA);
    }
    return dgGaussVectorBinary<OutType>(tA(), tB, eval, mutB);
}

struct DgGaussVecEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b) { return a + b; }
};

struct DgGaussVecEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b) { return a - b; }
};

struct DgGaussVecEvalCross
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a ^ b) { return a ^ b; }
};

struct DgGaussVecEvalDot
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b) { return a & b; }
};

struct DgGaussVecEvalOuter
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b) { return a * b; }
};

struct DgGaussVecEvalAndVT
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b) { return a & b; }
};

struct DgGaussVecEvalAndTV
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b) { return a & b; }
};

struct DgGaussVecEvalMag
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(mag(a)) { return mag(a); }
};

struct DgGaussVecEvalMagSqr
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(magSqr(a)) { return magSqr(a); }
};

struct DgGaussVecMutAddLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a + b; }
};

struct DgGaussVecMutAddRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a + b; }
};

struct DgGaussVecMutSubLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a - b; }
};

struct DgGaussVecMutSubRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

struct DgGaussVecMutCrossLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a ^ b; }
};

struct DgGaussVecMutCrossRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a ^ b; }
};

struct DgGaussVecMutAndLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a & b; }
};

struct DgGaussVecMutAndRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a & b; }
};

struct DgGaussVecMutNoop
{
    template<class A>
    inline void operator()(A&) const {}
};

struct DgGaussVecNoMutLeft
{
    template<class A, class B>
    inline void operator()(A&, const B&) const {}
};

struct DgGaussVecNoMutRight
{
    template<class A, class B>
    inline void operator()(const A&, B&) const {}
};

#define DG_GAUSS_VECTOR_DEFINE_UNARY(FUNC, OUTF, INF, EVAL, MUT)                            \
inline tmp<GaussField<OUTF>> FUNC (const GaussField<INF>& A)                                 \
{                                                                                              \
    return dgGaussVectorUnary<OUTF>(A, EVAL());                                               \
}                                                                                              \
inline tmp<GaussField<OUTF>> FUNC (const tmp<GaussField<INF>>& tA)                           \
{                                                                                              \
    return dgGaussVectorUnary<OUTF>(tA, EVAL(), MUT());                                       \
}

#define DG_GAUSS_VECTOR_DEFINE_BINARY(OP, OUTF, AF, BF, EVAL, MUTA, MUTB)                    \
inline tmp<GaussField<OUTF>> operator OP (const GaussField<AF>& A, const GaussField<BF>& B) \
{                                                                                              \
    return dgGaussVectorBinary<OUTF>(A, B, EVAL());                                           \
}                                                                                              \
inline tmp<GaussField<OUTF>> operator OP (const tmp<GaussField<AF>>& tA, const GaussField<BF>& B) \
{                                                                                              \
    return dgGaussVectorBinary<OUTF>(tA, B, EVAL(), MUTA());                                  \
}                                                                                              \
inline tmp<GaussField<OUTF>> operator OP (const GaussField<AF>& A, const tmp<GaussField<BF>>& tB) \
{                                                                                              \
    return dgGaussVectorBinary<OUTF>(A, tB, EVAL(), MUTB());                                  \
}                                                                                              \
inline tmp<GaussField<OUTF>> operator OP (const tmp<GaussField<AF>>& tA, const tmp<GaussField<BF>>& tB) \
{                                                                                              \
    return dgGaussVectorBinary<OUTF>(tA, tB, EVAL(), MUTA(), MUTB());                         \
}

// Vector - Vector
DG_GAUSS_VECTOR_DEFINE_BINARY(+, vector, vector, vector, DgGaussVecEvalAdd,   DgGaussVecMutAddLeft,   DgGaussVecMutAddRight)
DG_GAUSS_VECTOR_DEFINE_BINARY(-, vector, vector, vector, DgGaussVecEvalSub,   DgGaussVecMutSubLeft,   DgGaussVecMutSubRight)
DG_GAUSS_VECTOR_DEFINE_BINARY(^, vector, vector, vector, DgGaussVecEvalCross, DgGaussVecMutCrossLeft, DgGaussVecMutCrossRight)
DG_GAUSS_VECTOR_DEFINE_BINARY(&, scalar, vector, vector, DgGaussVecEvalDot,   DgGaussVecNoMutLeft,    DgGaussVecNoMutRight)
DG_GAUSS_VECTOR_DEFINE_BINARY(*, tensor, vector, vector, DgGaussVecEvalOuter, DgGaussVecNoMutLeft,    DgGaussVecNoMutRight)

DG_GAUSS_VECTOR_DEFINE_UNARY(mag,    scalar, vector, DgGaussVecEvalMag,    DgGaussVecMutNoop)
DG_GAUSS_VECTOR_DEFINE_UNARY(magSqr, scalar, vector, DgGaussVecEvalMagSqr, DgGaussVecMutNoop)

// Vector - Tensor
DG_GAUSS_VECTOR_DEFINE_BINARY(&, vector, vector, tensor, DgGaussVecEvalAndVT, DgGaussVecMutAndLeft, DgGaussVecNoMutRight)
DG_GAUSS_VECTOR_DEFINE_BINARY(&, vector, tensor, vector, DgGaussVecEvalAndTV, DgGaussVecNoMutLeft,  DgGaussVecMutAndRight)

// Vector - symmTensor
DG_GAUSS_VECTOR_DEFINE_BINARY(&, vector, vector, symmTensor, DgGaussVecEvalAndVT, DgGaussVecMutAndLeft, DgGaussVecNoMutRight)
DG_GAUSS_VECTOR_DEFINE_BINARY(&, vector, symmTensor, vector, DgGaussVecEvalAndTV, DgGaussVecNoMutLeft,  DgGaussVecMutAndRight)

// Vector - sphericalTensor
DG_GAUSS_VECTOR_DEFINE_BINARY(&, vector, vector, sphericalTensor, DgGaussVecEvalAndVT, DgGaussVecMutAndLeft, DgGaussVecNoMutRight)
DG_GAUSS_VECTOR_DEFINE_BINARY(&, vector, sphericalTensor, vector, DgGaussVecEvalAndTV, DgGaussVecNoMutLeft,  DgGaussVecMutAndRight)

#undef DG_GAUSS_VECTOR_DEFINE_UNARY
#undef DG_GAUSS_VECTOR_DEFINE_BINARY

} // End namespace Foam

#endif

// ************************************************************************* //
