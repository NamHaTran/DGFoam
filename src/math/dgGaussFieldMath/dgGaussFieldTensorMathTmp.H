/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Tensor math operations for GaussField (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldTensorMathTmp_H
#define Foam_dgGaussFieldTensorMathTmp_H

#include "GaussField.H"
#include "dgFaceGaussFieldTensorMathTmp.H"
#include "dgCellGaussFieldTensorMathTmp.H"
#include "Tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "scalar.H"
#include <type_traits>

namespace Foam
{

template<class OutType, class RefType>
inline tmp<GaussField<OutType>> dgGaussTensorNewLike(const GaussField<RefType>& A)
{
    return tmp<GaussField<OutType>>(new GaussField<OutType>(A.cellID(), A.dgMesh()));
}

template<class OutType, class InType, class Eval>
inline tmp<GaussField<OutType>> dgGaussTensorUnary
(
    const GaussField<InType>& A,
    const Eval& eval
)
{
    auto tRes = dgGaussTensorNewLike<OutType>(A);
    GaussField<OutType>& R = tRes.ref();
    R.cellField() = eval(A.cellField());
    R.faceField() = eval(A.faceField());
    return tRes;
}

template<class OutType, class InType, class Eval, class Mutate>
inline typename std::enable_if<std::is_same<OutType, InType>::value, tmp<GaussField<OutType>>>::type
dgGaussTensorUnary
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
    return dgGaussTensorUnary<OutType>(tA(), eval);
}

template<class OutType, class InType, class Eval, class Mutate>
inline typename std::enable_if<!std::is_same<OutType, InType>::value, tmp<GaussField<OutType>>>::type
dgGaussTensorUnary
(
    const tmp<GaussField<InType>>& tA,
    const Eval& eval,
    const Mutate&
)
{
    return dgGaussTensorUnary<OutType>(tA(), eval);
}

template<class OutType, class AType, class BType, class Eval>
inline tmp<GaussField<OutType>> dgGaussTensorBinary
(
    const GaussField<AType>& A,
    const GaussField<BType>& B,
    const Eval& eval
)
{
    auto tRes = dgGaussTensorNewLike<OutType>(A);
    GaussField<OutType>& R = tRes.ref();
    R.cellField() = eval(A.cellField(), B.cellField());
    R.faceField() = eval(A.faceField(), B.faceField());
    return tRes;
}

template<class OutType, class AType, class BType, class Eval, class MutA>
inline typename std::enable_if<std::is_same<OutType, AType>::value, tmp<GaussField<OutType>>>::type
dgGaussTensorBinary
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
    return dgGaussTensorBinary<OutType>(tA(), B, eval);
}

template<class OutType, class AType, class BType, class Eval, class MutA>
inline typename std::enable_if<!std::is_same<OutType, AType>::value, tmp<GaussField<OutType>>>::type
dgGaussTensorBinary
(
    const tmp<GaussField<AType>>& tA,
    const GaussField<BType>& B,
    const Eval& eval,
    const MutA&
)
{
    return dgGaussTensorBinary<OutType>(tA(), B, eval);
}

template<class OutType, class AType, class BType, class Eval, class MutB>
inline typename std::enable_if<std::is_same<OutType, BType>::value, tmp<GaussField<OutType>>>::type
dgGaussTensorBinary
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
    return dgGaussTensorBinary<OutType>(A, tB(), eval);
}

template<class OutType, class AType, class BType, class Eval, class MutB>
inline typename std::enable_if<!std::is_same<OutType, BType>::value, tmp<GaussField<OutType>>>::type
dgGaussTensorBinary
(
    const GaussField<AType>& A,
    const tmp<GaussField<BType>>& tB,
    const Eval& eval,
    const MutB&
)
{
    return dgGaussTensorBinary<OutType>(A, tB(), eval);
}

template<class OutType, class AType, class BType, class Eval, class MutA, class MutB>
inline tmp<GaussField<OutType>> dgGaussTensorBinary
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
        return dgGaussTensorBinary<OutType>(tA, tB(), eval, mutA);
    }
    return dgGaussTensorBinary<OutType>(tA(), tB, eval, mutB);
}

struct DgGaussTensorEvalDev
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(dev(a)) { return dev(a); }
};

struct DgGaussTensorEvalSymm
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(symm(a)) { return symm(a); }
};

struct DgGaussTensorEvalSkew
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(skew(a)) { return skew(a); }
};

struct DgGaussTensorEvalTr
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(tr(a)) { return tr(a); }
};

struct DgGaussTensorEvalDet
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(det(a)) { return det(a); }
};

struct DgGaussTensorEvalT
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(T(a)) { return T(a); }
};

struct DgGaussTensorEvalInv
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(inv(a)) { return inv(a); }
};

struct DgGaussTensorMutDev
{
    template<class A>
    inline void operator()(A& a) const { a = dev(a); }
};

struct DgGaussTensorMutSkew
{
    template<class A>
    inline void operator()(A& a) const { a = skew(a); }
};

struct DgGaussTensorMutT
{
    template<class A>
    inline void operator()(A& a) const { a = T(a); }
};

struct DgGaussTensorMutNoop
{
    template<class A>
    inline void operator()(A&) const {}
};

struct DgGaussTensorEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b) { return a + b; }
};

struct DgGaussTensorEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b) { return a - b; }
};

struct DgGaussTensorEvalAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b) { return a & b; }
};

struct DgGaussTensorEvalAndRev
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(b & a) { return b & a; }
};

struct DgGaussTensorEvalAndAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a && b) { return a && b; }
};

struct DgGaussTensorEvalAndAndRev
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(b && a) { return b && a; }
};

struct DgGaussTensorMutAddLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a + b; }
};

struct DgGaussTensorMutAddRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a + b; }
};

struct DgGaussTensorMutSubLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a - b; }
};

struct DgGaussTensorMutSubRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

struct DgGaussTensorMutAndLeft
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a & b; }
};

struct DgGaussTensorMutAndRight
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a & b; }
};

struct DgGaussTensorMutAndRightRev
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = b & a; }
};

struct DgGaussTensorNoMutLeft
{
    template<class A, class B>
    inline void operator()(A&, const B&) const {}
};

struct DgGaussTensorNoMutRight
{
    template<class A, class B>
    inline void operator()(const A&, B&) const {}
};

#define DG_GAUSS_TENSOR_DEFINE_UNARY(FUNC, OUTF, INF, EVAL, MUT)                              \
inline tmp<GaussField<OUTF>> FUNC (const GaussField<INF>& A)                                  \
{                                                                                               \
    return dgGaussTensorUnary<OUTF>(A, EVAL());                                                \
}                                                                                               \
inline tmp<GaussField<OUTF>> FUNC (const tmp<GaussField<INF>>& tA)                            \
{                                                                                               \
    return dgGaussTensorUnary<OUTF>(tA, EVAL(), MUT());                                        \
}

#define DG_GAUSS_TENSOR_DEFINE_BINARY(OP, OUTF, AF, BF, EVAL, MUTA, MUTB)                     \
inline tmp<GaussField<OUTF>> operator OP (const GaussField<AF>& A, const GaussField<BF>& B)  \
{                                                                                               \
    return dgGaussTensorBinary<OUTF>(A, B, EVAL());                                            \
}                                                                                               \
inline tmp<GaussField<OUTF>> operator OP (const tmp<GaussField<AF>>& tA, const GaussField<BF>& B) \
{                                                                                               \
    return dgGaussTensorBinary<OUTF>(tA, B, EVAL(), MUTA());                                   \
}                                                                                               \
inline tmp<GaussField<OUTF>> operator OP (const GaussField<AF>& A, const tmp<GaussField<BF>>& tB) \
{                                                                                               \
    return dgGaussTensorBinary<OUTF>(A, tB, EVAL(), MUTB());                                   \
}                                                                                               \
inline tmp<GaussField<OUTF>> operator OP (const tmp<GaussField<AF>>& tA, const tmp<GaussField<BF>>& tB) \
{                                                                                               \
    return dgGaussTensorBinary<OUTF>(tA, tB, EVAL(), MUTA(), MUTB());                          \
}

// tensor utilities
DG_GAUSS_TENSOR_DEFINE_UNARY(dev,  tensor,     tensor,     DgGaussTensorEvalDev,   DgGaussTensorMutDev)
DG_GAUSS_TENSOR_DEFINE_UNARY(symm, symmTensor, tensor,     DgGaussTensorEvalSymm,  DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(skew, tensor,     tensor,     DgGaussTensorEvalSkew,  DgGaussTensorMutSkew)
DG_GAUSS_TENSOR_DEFINE_UNARY(T,    tensor,     tensor,     DgGaussTensorEvalT,     DgGaussTensorMutT)
DG_GAUSS_TENSOR_DEFINE_UNARY(tr,   scalar,     tensor,     DgGaussTensorEvalTr,    DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(det,  scalar,     tensor,     DgGaussTensorEvalDet,   DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(inv,  tensor,     tensor,     DgGaussTensorEvalInv,   DgGaussTensorMutNoop)

// symmTensor utilities
DG_GAUSS_TENSOR_DEFINE_UNARY(tr,   scalar,     symmTensor, DgGaussTensorEvalTr,    DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(dev,  symmTensor, symmTensor, DgGaussTensorEvalDev,   DgGaussTensorMutDev)
DG_GAUSS_TENSOR_DEFINE_UNARY(inv,  tensor,     symmTensor, DgGaussTensorEvalInv,   DgGaussTensorMutNoop)

// tensor - tensor
DG_GAUSS_TENSOR_DEFINE_BINARY(+,  tensor, tensor,     tensor,     DgGaussTensorEvalAdd,    DgGaussTensorMutAddLeft, DgGaussTensorMutAddRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(-,  tensor, tensor,     tensor,     DgGaussTensorEvalSub,    DgGaussTensorMutSubLeft, DgGaussTensorMutSubRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&,  tensor, tensor,     tensor,     DgGaussTensorEvalAnd,    DgGaussTensorMutAndLeft, DgGaussTensorMutAndRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&&, scalar, tensor,     tensor,     DgGaussTensorEvalAndAnd, DgGaussTensorNoMutLeft,  DgGaussTensorNoMutRight)

// tensor - symmTensor
DG_GAUSS_TENSOR_DEFINE_BINARY(+,  tensor, tensor,     symmTensor, DgGaussTensorEvalAdd,    DgGaussTensorMutAddLeft, DgGaussTensorNoMutRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(-,  tensor, tensor,     symmTensor, DgGaussTensorEvalSub,    DgGaussTensorMutSubLeft, DgGaussTensorNoMutRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&,  tensor, tensor,     symmTensor, DgGaussTensorEvalAnd,    DgGaussTensorMutAndLeft, DgGaussTensorNoMutRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&&, scalar, tensor,     symmTensor, DgGaussTensorEvalAndAnd, DgGaussTensorNoMutLeft,  DgGaussTensorNoMutRight)

// symmTensor - tensor (keep original semantics for &, && as B & A, B && A)
DG_GAUSS_TENSOR_DEFINE_BINARY(+,  tensor, symmTensor, tensor,     DgGaussTensorEvalAdd,    DgGaussTensorNoMutLeft,  DgGaussTensorMutAddRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(-,  tensor, symmTensor, tensor,     DgGaussTensorEvalSub,    DgGaussTensorNoMutLeft,  DgGaussTensorMutSubRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&,  tensor, symmTensor, tensor,     DgGaussTensorEvalAndRev, DgGaussTensorNoMutLeft,  DgGaussTensorMutAndRightRev)
DG_GAUSS_TENSOR_DEFINE_BINARY(&&, scalar, symmTensor, tensor,     DgGaussTensorEvalAndAndRev, DgGaussTensorNoMutLeft, DgGaussTensorNoMutRight)

// symmTensor - symmTensor
DG_GAUSS_TENSOR_DEFINE_BINARY(+,  symmTensor, symmTensor, symmTensor, DgGaussTensorEvalAdd,    DgGaussTensorMutAddLeft, DgGaussTensorMutAddRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(-,  symmTensor, symmTensor, symmTensor, DgGaussTensorEvalSub,    DgGaussTensorMutSubLeft, DgGaussTensorMutSubRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&,  tensor,     symmTensor, symmTensor, DgGaussTensorEvalAnd,    DgGaussTensorNoMutLeft,  DgGaussTensorNoMutRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&&, scalar,     symmTensor, symmTensor, DgGaussTensorEvalAndAnd, DgGaussTensorNoMutLeft,  DgGaussTensorNoMutRight)

#undef DG_GAUSS_TENSOR_DEFINE_UNARY
#undef DG_GAUSS_TENSOR_DEFINE_BINARY

} // End namespace Foam

#endif

// ************************************************************************* //
