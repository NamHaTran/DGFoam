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

/** \brief Allocate output GaussField with geometry copied from reference field. */
template<class OutType, class RefType>
inline tmp<GaussField<OutType>> dgGaussTensorNewLike(const GaussField<RefType>& A)
{
    return tmp<GaussField<OutType>>(new GaussField<OutType>(A.cellID(), A.dgMesh()));
}

/** \brief Core unary dispatcher for non-tmp tensor GaussField input. */
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

/** \brief Core binary dispatcher for non-tmp tensor GaussField inputs. */
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

/** \brief Evaluate `dev(.)` at one Gauss point. */
struct DgGaussTensorEvalDev
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(dev(a)) { return dev(a); }
};

/** \brief Evaluate `symm(.)` at one Gauss point. */
struct DgGaussTensorEvalSymm
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(symm(a)) { return symm(a); }
};

/** \brief Evaluate `skew(.)` at one Gauss point. */
struct DgGaussTensorEvalSkew
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(skew(a)) { return skew(a); }
};

/** \brief Evaluate `tr(.)` at one Gauss point. */
struct DgGaussTensorEvalTr
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(tr(a)) { return tr(a); }
};

/** \brief Evaluate `det(.)` at one Gauss point. */
struct DgGaussTensorEvalDet
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(det(a)) { return det(a); }
};

/** \brief Evaluate transpose `T(.)` at one Gauss point. */
struct DgGaussTensorEvalT
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(T(a)) { return T(a); }
};

/** \brief Evaluate `inv(.)` at one Gauss point. */
struct DgGaussTensorEvalInv
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(inv(a)) { return inv(a); }
};

/** \brief In-place mutate tensor by `dev(.)`. */
struct DgGaussTensorMutDev
{
    template<class A>
    inline void operator()(A& a) const { a = dev(a); }
};

/** \brief In-place mutate tensor by `skew(.)`. */
struct DgGaussTensorMutSkew
{
    template<class A>
    inline void operator()(A& a) const { a = skew(a); }
};

/** \brief In-place mutate tensor by transpose `T(.)`. */
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

/** \brief Evaluate binary addition at one Gauss point. */
struct DgGaussTensorEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b) { return a + b; }
};

/** \brief Evaluate binary subtraction at one Gauss point. */
struct DgGaussTensorEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b) { return a - b; }
};

/** \brief Evaluate single contraction `&` at one Gauss point. */
struct DgGaussTensorEvalAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b) { return a & b; }
};

/** \brief Evaluate reverse contraction `b & a` at one Gauss point. */
struct DgGaussTensorEvalAndRev
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(b & a) { return b & a; }
};

/** \brief Evaluate double contraction `&&` at one Gauss point. */
struct DgGaussTensorEvalAndAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a && b) { return a && b; }
};

/** \brief Evaluate reverse double contraction `b && a` at one Gauss point. */
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

/** \brief Unary tensor utilities: `dev`, `symm`, `skew`, `T`, `tr`, `det`, `inv`. */
// tensor utilities
DG_GAUSS_TENSOR_DEFINE_UNARY(dev,  tensor,     tensor,     DgGaussTensorEvalDev,   DgGaussTensorMutDev)
DG_GAUSS_TENSOR_DEFINE_UNARY(symm, symmTensor, tensor,     DgGaussTensorEvalSymm,  DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(skew, tensor,     tensor,     DgGaussTensorEvalSkew,  DgGaussTensorMutSkew)
DG_GAUSS_TENSOR_DEFINE_UNARY(T,    tensor,     tensor,     DgGaussTensorEvalT,     DgGaussTensorMutT)
DG_GAUSS_TENSOR_DEFINE_UNARY(tr,   scalar,     tensor,     DgGaussTensorEvalTr,    DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(det,  scalar,     tensor,     DgGaussTensorEvalDet,   DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(inv,  tensor,     tensor,     DgGaussTensorEvalInv,   DgGaussTensorMutNoop)

/** \brief Unary symmTensor utilities: `tr`, `dev`, `inv`. */
// symmTensor utilities
DG_GAUSS_TENSOR_DEFINE_UNARY(tr,   scalar,     symmTensor, DgGaussTensorEvalTr,    DgGaussTensorMutNoop)
DG_GAUSS_TENSOR_DEFINE_UNARY(dev,  symmTensor, symmTensor, DgGaussTensorEvalDev,   DgGaussTensorMutDev)
DG_GAUSS_TENSOR_DEFINE_UNARY(inv,  tensor,     symmTensor, DgGaussTensorEvalInv,   DgGaussTensorMutNoop)

/** \brief Binary tensor-tensor operators: `+`, `-`, `&`, `&&`. */
// tensor - tensor
DG_GAUSS_TENSOR_DEFINE_BINARY(+,  tensor, tensor,     tensor,     DgGaussTensorEvalAdd,    DgGaussTensorMutAddLeft, DgGaussTensorMutAddRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(-,  tensor, tensor,     tensor,     DgGaussTensorEvalSub,    DgGaussTensorMutSubLeft, DgGaussTensorMutSubRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&,  tensor, tensor,     tensor,     DgGaussTensorEvalAnd,    DgGaussTensorMutAndLeft, DgGaussTensorMutAndRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&&, scalar, tensor,     tensor,     DgGaussTensorEvalAndAnd, DgGaussTensorNoMutLeft,  DgGaussTensorNoMutRight)

/** \brief Binary tensor-symmTensor operators: `+`, `-`, `&`, `&&`. */
// tensor - symmTensor
DG_GAUSS_TENSOR_DEFINE_BINARY(+,  tensor, tensor,     symmTensor, DgGaussTensorEvalAdd,    DgGaussTensorMutAddLeft, DgGaussTensorNoMutRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(-,  tensor, tensor,     symmTensor, DgGaussTensorEvalSub,    DgGaussTensorMutSubLeft, DgGaussTensorNoMutRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&,  tensor, tensor,     symmTensor, DgGaussTensorEvalAnd,    DgGaussTensorMutAndLeft, DgGaussTensorNoMutRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&&, scalar, tensor,     symmTensor, DgGaussTensorEvalAndAnd, DgGaussTensorNoMutLeft,  DgGaussTensorNoMutRight)

/** \brief Binary symmTensor-tensor operators: `+`, `-`, `&`, `&&`. */
// symmTensor - tensor (keep original semantics for &, && as B & A, B && A)
DG_GAUSS_TENSOR_DEFINE_BINARY(+,  tensor, symmTensor, tensor,     DgGaussTensorEvalAdd,    DgGaussTensorNoMutLeft,  DgGaussTensorMutAddRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(-,  tensor, symmTensor, tensor,     DgGaussTensorEvalSub,    DgGaussTensorNoMutLeft,  DgGaussTensorMutSubRight)
DG_GAUSS_TENSOR_DEFINE_BINARY(&,  tensor, symmTensor, tensor,     DgGaussTensorEvalAndRev, DgGaussTensorNoMutLeft,  DgGaussTensorMutAndRightRev)
DG_GAUSS_TENSOR_DEFINE_BINARY(&&, scalar, symmTensor, tensor,     DgGaussTensorEvalAndAndRev, DgGaussTensorNoMutLeft, DgGaussTensorNoMutRight)

/** \brief Binary symmTensor-symmTensor operators: `+`, `-`, `&`, `&&`. */
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
