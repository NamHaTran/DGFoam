/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
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

// No-op mutators for operations that cannot/should not reuse tmp storage
struct dgNoMutateLeft
{
    template<class LhsElem, class RhsElem>
    inline void operator()(LhsElem&, const RhsElem&) const {}
};

struct dgNoMutateRight
{
    template<class LhsElem, class RhsElem>
    inline void operator()(const LhsElem&, RhsElem&) const {}
};

struct dgNoMutateUnary
{
    template<class Elem>
    inline void operator()(Elem&) const {}
};

template<class FieldType>
inline tmp<FieldType> dgNewLike(const FieldType& A)
{
    return tmp<FieldType>(new FieldType(A.cellID(), A.dgMesh()));
}

// --------------------------- Unary core --------------------------- //

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

// --------------------------- Binary core -------------------------- //

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

// ------------------------- Compare core --------------------------- //

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


// ---------------------------- Macros ------------------------------ //
// Eval:    ElemOut operator()(const ElemA&, const ElemB&) const;
// MutateA: void operator()(ElemA&, const ElemB&) const;
// MutateB: void operator()(const ElemA&, ElemB&) const;

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

// Eval: ElemOut operator()(const ElemIn&) const;
// Mut:  void operator()(ElemIn&) const;
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

