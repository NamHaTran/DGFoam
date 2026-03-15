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
    Core helpers/macros to reduce overload boilerplate for tmp/non-tmp
    arithmetic on index-based Gauss fields (eg, cellGaussField).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldTmpOpCore_H
#define Foam_dgGaussFieldTmpOpCore_H

#include "tmp.H"
#include <type_traits>

namespace Foam
{

/**
 * \brief Quick guide to add a new tmp-aware operator using this core.
 *
 * Steps:
 * 1. Define `Eval` functor for element-wise computation.
 * 2. Define optional `Mutate` functors for in-place tmp reuse.
 * 3. Use macro generators at file end:
 *    - `DG_TMP_DEFINE_BINARY_OPERATOR`
 *    - `DG_TMP_DEFINE_UNARY_FUNCTION`
 *    - `DG_TMP_DEFINE_COMPARE_OPERATOR`
 * 4. For mixed scalar-field operators, reuse helpers:
 *    - `dgTmpBinary`, `dgTmpUnary`, `dgTmpCompare`
 *    - or the field-scalar helpers in specialized cores.
 * 5. Use `dgNoMutateLeft/Right/Unary` when in-place reuse is not valid.
 *
 * Core function roles:
 * - `dgNewLike`: allocate result field with same geometry metadata.
 * - `dgTmpUnary`: unary kernels for non-tmp/tmp inputs.
 * - `dgTmpBinary`: binary kernels for all tmp/non-tmp combinations.
 * - `dgTmpCompare`: element-wise predicate checks with early exit.
 */

/** \brief No-op mutator for lhs when tmp reuse is disabled. */
struct dgNoMutateLeft
{
    template<class LhsElem, class RhsElem>
    inline void operator()(LhsElem&, const RhsElem&) const {}
};

/** \brief No-op mutator for rhs when tmp reuse is disabled. */
struct dgNoMutateRight
{
    template<class LhsElem, class RhsElem>
    inline void operator()(const LhsElem&, RhsElem&) const {}
};

/** \brief No-op mutator for unary tmp reuse-disabled paths. */
struct dgNoMutateUnary
{
    template<class Elem>
    inline void operator()(Elem&) const {}
};

/** \brief Allocate output field with mesh/cell metadata copied from reference field. */
template<class FieldType>
inline tmp<FieldType> dgNewLike(const FieldType& A)
{
    return tmp<FieldType>(new FieldType(A.cellID(), A.dgMesh()));
}

/** \brief Unary core for non-tmp input fields. */

template<class OutField, class InField, class Eval>
inline tmp<OutField> dgTmpUnary
(
    const InField& A,
    const Eval& eval
)
{
    auto tRes = dgNewLike<OutField>(A);
    OutField& R = tRes.ref();

    for (label i = 0; i < R.size(); ++i)
    {
        R[i] = eval(A[i]);
    }

    return tRes;
}

/** \brief Unary core for tmp inputs with possible in-place reuse. */
template<class OutField, class InField, class Eval, class Mutate>
inline typename std::enable_if<std::is_same<OutField, InField>::value, tmp<OutField>>::type
dgTmpUnary
(
    const tmp<InField>& tA,
    const Eval& eval,
    const Mutate& mutate
)
{
    if (tA.isTmp())
    {
        InField& A = const_cast<tmp<InField>&>(tA).ref();

        for (label i = 0; i < A.size(); ++i)
        {
            mutate(A[i]);
        }

        return tA;
    }

    return dgTmpUnary<OutField>(tA(), eval);
}

/** \brief Unary core for tmp inputs when in-place reuse is not possible. */
template<class OutField, class InField, class Eval, class Mutate>
inline typename std::enable_if<!std::is_same<OutField, InField>::value, tmp<OutField>>::type
dgTmpUnary
(
    const tmp<InField>& tA,
    const Eval& eval,
    const Mutate&
)
{
    return dgTmpUnary<OutField>(tA(), eval);
}

/** \brief Binary core for non-tmp/non-tmp inputs. */

template<class OutField, class AField, class BField, class Eval>
inline tmp<OutField> dgTmpBinary
(
    const AField& A,
    const BField& B,
    const Eval& eval
)
{
    auto tRes = dgNewLike<OutField>(A);
    OutField& R = tRes.ref();

    for (label i = 0; i < R.size(); ++i)
    {
        R[i] = eval(A[i], B[i]);
    }

    return tRes;
}

/** \brief Binary core for tmp lhs + non-tmp rhs with optional lhs reuse. */
template<class OutField, class AField, class BField, class Eval, class MutateA>
inline typename std::enable_if<std::is_same<OutField, AField>::value, tmp<OutField>>::type
dgTmpBinary
(
    const tmp<AField>& tA,
    const BField& B,
    const Eval& eval,
    const MutateA& mutateA
)
{
    if (tA.isTmp())
    {
        AField& A = const_cast<tmp<AField>&>(tA).ref();

        for (label i = 0; i < A.size(); ++i)
        {
            mutateA(A[i], B[i]);
        }

        return tA;
    }

    return dgTmpBinary<OutField>(tA(), B, eval);
}

/** \brief Binary core for tmp lhs + non-tmp rhs without lhs reuse. */
template<class OutField, class AField, class BField, class Eval, class MutateA>
inline typename std::enable_if<!std::is_same<OutField, AField>::value, tmp<OutField>>::type
dgTmpBinary
(
    const tmp<AField>& tA,
    const BField& B,
    const Eval& eval,
    const MutateA&
)
{
    return dgTmpBinary<OutField>(tA(), B, eval);
}

/** \brief Binary core for non-tmp lhs + tmp rhs with optional rhs reuse. */
template<class OutField, class AField, class BField, class Eval, class MutateB>
inline typename std::enable_if<std::is_same<OutField, BField>::value, tmp<OutField>>::type
dgTmpBinary
(
    const AField& A,
    const tmp<BField>& tB,
    const Eval& eval,
    const MutateB& mutateB
)
{
    if (tB.isTmp())
    {
        BField& B = const_cast<tmp<BField>&>(tB).ref();

        for (label i = 0; i < B.size(); ++i)
        {
            mutateB(A[i], B[i]);
        }

        return tB;
    }

    return dgTmpBinary<OutField>(A, tB(), eval);
}

/** \brief Binary core for non-tmp lhs + tmp rhs without rhs reuse. */
template<class OutField, class AField, class BField, class Eval, class MutateB>
inline typename std::enable_if<!std::is_same<OutField, BField>::value, tmp<OutField>>::type
dgTmpBinary
(
    const AField& A,
    const tmp<BField>& tB,
    const Eval& eval,
    const MutateB&
)
{
    return dgTmpBinary<OutField>(A, tB(), eval);
}

/** \brief Binary core for tmp/tmp inputs, preferring lhs reuse path first. */
template<class OutField, class AField, class BField, class Eval, class MutateA, class MutateB>
inline tmp<OutField> dgTmpBinary
(
    const tmp<AField>& tA,
    const tmp<BField>& tB,
    const Eval& eval,
    const MutateA& mutateA,
    const MutateB& mutateB
)
{
    if (tA.isTmp())
    {
        return dgTmpBinary<OutField>(tA, tB(), eval, mutateA);
    }

    return dgTmpBinary<OutField>(tA(), tB, eval, mutateB);
}

/** \brief Compare two non-tmp fields element-wise with early exit. */

/** \brief Compare tmp lhs against non-tmp rhs. */
template<class AField, class BField, class Pred>
inline bool dgTmpCompare
(
    const AField& A,
    const BField& B,
    const Pred& pred
)
{
    for (label i = 0; i < A.size(); ++i)
    {
        if (!pred(A[i], B[i]))
        {
            return false;
        }
    }

    return true;
}

/** \brief Compare non-tmp lhs against tmp rhs. */
template<class AField, class BField, class Pred>
inline bool dgTmpCompare
(
    const tmp<AField>& tA,
    const BField& B,
    const Pred& pred
)
{
    return dgTmpCompare(tA(), B, pred);
}

/** \brief Compare tmp lhs against tmp rhs. */
template<class AField, class BField, class Pred>
inline bool dgTmpCompare
(
    const AField& A,
    const tmp<BField>& tB,
    const Pred& pred
)
{
    return dgTmpCompare(A, tB(), pred);
}

template<class AField, class BField, class Pred>
inline bool dgTmpCompare
(
    const tmp<AField>& tA,
    const tmp<BField>& tB,
    const Pred& pred
)
{
    return dgTmpCompare(tA(), tB(), pred);
}

} // End namespace Foam


/** \brief Macro family to generate binary operator overloads for all tmp/non-tmp combinations. */

/** \brief Generate binary operator overloads for all tmp/non-tmp combinations. */
#define DG_TMP_DEFINE_BINARY_OPERATOR(OP, OUTFIELD, AFIELD, BFIELD, EVAL, MUTA, MUTB) \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const BFIELD& B)                     \
{                                                                                         \
    return dgTmpBinary<OUTFIELD>(A, B, EVAL());                                           \
}                                                                                         \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const BFIELD& B)               \
{                                                                                         \
    return dgTmpBinary<OUTFIELD>(tA, B, EVAL(), MUTA());                                  \
}                                                                                         \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const tmp<BFIELD>& tB)               \
{                                                                                         \
    return dgTmpBinary<OUTFIELD>(A, tB, EVAL(), MUTB());                                  \
}                                                                                         \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const tmp<BFIELD>& tB)         \
{                                                                                         \
    return dgTmpBinary<OUTFIELD>(tA, tB, EVAL(), MUTA(), MUTB());                         \
}

/** \brief Macro family to generate unary function overloads for non-tmp/tmp inputs. */
/** \brief Generate unary function overloads for non-tmp/tmp inputs. */
#define DG_TMP_DEFINE_UNARY_FUNCTION(FUNC, OUTFIELD, INFIELD, EVAL, MUT) \
inline tmp<OUTFIELD> FUNC (const INFIELD& A)                              \
{                                                                          \
    return dgTmpUnary<OUTFIELD>(A, EVAL());                                \
}                                                                          \
inline tmp<OUTFIELD> FUNC (const tmp<INFIELD>& tA)                        \
{                                                                          \
    return dgTmpUnary<OUTFIELD>(tA, EVAL(), MUT());                        \
}

// Pred: bool operator()(const ElemA&, const ElemB&) const;
/** \brief Generate `(field OP field)` comparison overloads with tmp-aware path. */
#define DG_TMP_DEFINE_COMPARE_OPERATOR(OP, AFIELD, BFIELD, PRED) \
inline bool operator OP (const AFIELD& A, const BFIELD& B)       \
{                                                                 \
    return dgTmpCompare(A, B, PRED());                            \
}                                                                 \
inline bool operator OP (const tmp<AFIELD>& tA, const BFIELD& B) \
{                                                                 \
    return dgTmpCompare(tA, B, PRED());                           \
}                                                                 \
inline bool operator OP (const AFIELD& A, const tmp<BFIELD>& tB) \
{                                                                 \
    return dgTmpCompare(A, tB, PRED());                           \
}                                                                 \
inline bool operator OP (const tmp<AFIELD>& tA, const tmp<BFIELD>& tB) \
{                                                                      \
    return dgTmpCompare(tA, tB, PRED());                               \
}

#endif
