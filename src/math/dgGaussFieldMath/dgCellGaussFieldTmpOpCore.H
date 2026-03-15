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
    Core helpers/macros to reduce tmp/non-tmp overload boilerplate for
    cellGaussField-style indexed fields.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldTmpOpCore_H
#define Foam_dgCellGaussFieldTmpOpCore_H

#include "tmp.H"
#include <type_traits>

namespace Foam
{
namespace dgCellTmp
{

/**
 * \brief Quick guide to add a new operator using this core.
 *
 * Steps:
 * 1. Define an `Eval` functor that computes element-wise result.
 * 2. Define optional `Mutate` functors for in-place reuse when output type
 *    matches lhs/rhs tmp storage.
 * 3. Choose proper generator macro at file end:
 *    - `DG_CELL_TMP_DEFINE_BINARY_OPERATOR`
 *    - `DG_CELL_TMP_DEFINE_UNARY_FUNCTION`
 *    - `DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR`
 *    - `DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR`
 *    - `DG_CELL_TMP_DEFINE_FIELD_SCALAR_FUNCTION`
 *    - comparison macros for boolean operators
 * 4. Use `NoMutateLeft/Right/Unary` for non-reusable paths.
 *
 * Core function roles:
 * - `newLike`: allocate result field with same mesh/cell metadata.
 * - `unary`: apply unary eval for non-tmp/tmp inputs.
 * - `binary`: apply binary eval for all tmp/non-tmp combinations.
 * - `binaryFieldScalar`/`binaryScalarField`: mixed field-scalar kernels.
 * - `compare`/`compareFieldScalar`/`compareScalarField`: element-wise predicates.
 */

struct NoMutateLeft
{
    template<class LhsElem, class RhsElem>
    inline void operator()(LhsElem&, const RhsElem&) const {}
};

struct NoMutateRight
{
    template<class LhsElem, class RhsElem>
    inline void operator()(const LhsElem&, RhsElem&) const {}
};

struct NoMutateUnary
{
    template<class Elem>
    inline void operator()(Elem&) const {}
};

/** \brief Allocate output field with geometry copied from reference field. */
template<class OutField, class RefField>
inline tmp<OutField> newLike(const RefField& A)
{
    return tmp<OutField>(new OutField(A.cellID(), A.dgMesh()));
}

/** \brief Unary kernel for non-tmp input field. */
template<class OutField, class InField, class Eval>
inline tmp<OutField> unary
(
    const InField& A,
    const Eval& eval
)
{
    auto tRes = newLike<OutField>(A);
    OutField& R = tRes.ref();

    for (label i = 0; i < R.size(); ++i)
    {
        R[i] = eval(A[i]);
    }

    return tRes;
}

/** \brief Unary kernel for tmp input with in-place reuse when possible. */
template<class OutField, class InField, class Eval, class Mutate>
inline typename std::enable_if<std::is_same<OutField, InField>::value, tmp<OutField>>::type
unary
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

    return unary<OutField>(tA(), eval);
}

/** \brief Unary kernel for tmp input when in-place reuse is not possible. */
template<class OutField, class InField, class Eval, class Mutate>
inline typename std::enable_if<!std::is_same<OutField, InField>::value, tmp<OutField>>::type
unary
(
    const tmp<InField>& tA,
    const Eval& eval,
    const Mutate&
)
{
    return unary<OutField>(tA(), eval);
}

/** \brief Binary kernel for non-tmp/non-tmp inputs. */
template<class OutField, class AField, class BField, class Eval>
inline tmp<OutField> binary
(
    const AField& A,
    const BField& B,
    const Eval& eval
)
{
    auto tRes = newLike<OutField>(A);
    OutField& R = tRes.ref();

    for (label i = 0; i < R.size(); ++i)
    {
        R[i] = eval(A[i], B[i]);
    }

    return tRes;
}

/** \brief Binary kernel for tmp lhs + non-tmp rhs with optional lhs reuse. */
template<class OutField, class AField, class BField, class Eval, class MutateA>
inline typename std::enable_if<std::is_same<OutField, AField>::value, tmp<OutField>>::type
binary
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

    return binary<OutField>(tA(), B, eval);
}

template<class OutField, class AField, class BField, class Eval, class MutateA>
inline typename std::enable_if<!std::is_same<OutField, AField>::value, tmp<OutField>>::type
binary
(
    const tmp<AField>& tA,
    const BField& B,
    const Eval& eval,
    const MutateA&
)
{
    return binary<OutField>(tA(), B, eval);
}

/** \brief Binary kernel for non-tmp lhs + tmp rhs with optional rhs reuse. */
template<class OutField, class AField, class BField, class Eval, class MutateB>
inline typename std::enable_if<std::is_same<OutField, BField>::value, tmp<OutField>>::type
binary
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

    return binary<OutField>(A, tB(), eval);
}

template<class OutField, class AField, class BField, class Eval, class MutateB>
inline typename std::enable_if<!std::is_same<OutField, BField>::value, tmp<OutField>>::type
binary
(
    const AField& A,
    const tmp<BField>& tB,
    const Eval& eval,
    const MutateB&
)
{
    return binary<OutField>(A, tB(), eval);
}

template<class OutField, class AField, class BField, class Eval, class MutateA, class MutateB>
inline tmp<OutField> binary
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
        return binary<OutField>(tA, tB(), eval, mutateA);
    }

    return binary<OutField>(tA(), tB, eval, mutateB);
}

/** \brief Field-scalar kernel for non-tmp field input. */
template<class OutField, class AField, class ScalarType, class Eval>
inline tmp<OutField> binaryFieldScalar
(
    const AField& A,
    const ScalarType& b,
    const Eval& eval
)
{
    auto tRes = newLike<OutField>(A);
    OutField& R = tRes.ref();

    for (label i = 0; i < R.size(); ++i)
    {
        R[i] = eval(A[i], b);
    }

    return tRes;
}

template<class OutField, class AField, class ScalarType, class Eval, class MutateA>
inline typename std::enable_if<std::is_same<OutField, AField>::value, tmp<OutField>>::type
binaryFieldScalar
(
    const tmp<AField>& tA,
    const ScalarType& b,
    const Eval& eval,
    const MutateA& mutateA
)
{
    if (tA.isTmp())
    {
        AField& A = const_cast<tmp<AField>&>(tA).ref();

        for (label i = 0; i < A.size(); ++i)
        {
            mutateA(A[i], b);
        }

        return tA;
    }

    return binaryFieldScalar<OutField>(tA(), b, eval);
}

template<class OutField, class AField, class ScalarType, class Eval, class MutateA>
inline typename std::enable_if<!std::is_same<OutField, AField>::value, tmp<OutField>>::type
binaryFieldScalar
(
    const tmp<AField>& tA,
    const ScalarType& b,
    const Eval& eval,
    const MutateA&
)
{
    return binaryFieldScalar<OutField>(tA(), b, eval);
}

/** \brief Scalar-field kernel for non-tmp field input. */
template<class OutField, class ScalarType, class BField, class Eval>
inline tmp<OutField> binaryScalarField
(
    const ScalarType& a,
    const BField& B,
    const Eval& eval
)
{
    auto tRes = newLike<OutField>(B);
    OutField& R = tRes.ref();

    for (label i = 0; i < R.size(); ++i)
    {
        R[i] = eval(a, B[i]);
    }

    return tRes;
}

template<class OutField, class ScalarType, class BField, class Eval, class MutateB>
inline typename std::enable_if<std::is_same<OutField, BField>::value, tmp<OutField>>::type
binaryScalarField
(
    const ScalarType& a,
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
            mutateB(a, B[i]);
        }

        return tB;
    }

    return binaryScalarField<OutField>(a, tB(), eval);
}

template<class OutField, class ScalarType, class BField, class Eval, class MutateB>
inline typename std::enable_if<!std::is_same<OutField, BField>::value, tmp<OutField>>::type
binaryScalarField
(
    const ScalarType& a,
    const tmp<BField>& tB,
    const Eval& eval,
    const MutateB&
)
{
    return binaryScalarField<OutField>(a, tB(), eval);
}

/** \brief Compare two fields element-wise with early exit. */
template<class AField, class BField, class Pred>
inline bool compare
(
    const AField& A,
    const BField& B,
    const Pred& pred
)
{
    if (A.size() != B.size())
    {
        return false;
    }

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
inline bool compare
(
    const tmp<AField>& tA,
    const BField& B,
    const Pred& pred
)
{
    return compare(tA(), B, pred);
}

template<class AField, class BField, class Pred>
inline bool compare
(
    const AField& A,
    const tmp<BField>& tB,
    const Pred& pred
)
{
    return compare(A, tB(), pred);
}

template<class AField, class BField, class Pred>
inline bool compare
(
    const tmp<AField>& tA,
    const tmp<BField>& tB,
    const Pred& pred
)
{
    return compare(tA(), tB(), pred);
}

/** \brief Compare each field value against one scalar. */
template<class AField, class ScalarType, class Pred>
inline bool compareFieldScalar
(
    const AField& A,
    const ScalarType& b,
    const Pred& pred
)
{
    for (label i = 0; i < A.size(); ++i)
    {
        if (!pred(A[i], b))
        {
            return false;
        }
    }

    return true;
}

template<class AField, class ScalarType, class Pred>
inline bool compareFieldScalar
(
    const tmp<AField>& tA,
    const ScalarType& b,
    const Pred& pred
)
{
    return compareFieldScalar(tA(), b, pred);
}

/** \brief Compare one scalar against each field value. */
template<class ScalarType, class BField, class Pred>
inline bool compareScalarField
(
    const ScalarType& a,
    const BField& B,
    const Pred& pred
)
{
    for (label i = 0; i < B.size(); ++i)
    {
        if (!pred(a, B[i]))
        {
            return false;
        }
    }

    return true;
}

template<class ScalarType, class BField, class Pred>
inline bool compareScalarField
(
    const ScalarType& a,
    const tmp<BField>& tB,
    const Pred& pred
)
{
    return compareScalarField(a, tB(), pred);
}

} // End namespace dgCellTmp
} // End namespace Foam

/** \brief Generate binary operator overloads for all tmp/non-tmp combinations. */
#define DG_CELL_TMP_DEFINE_BINARY_OPERATOR(OP, OUTFIELD, AFIELD, BFIELD, EVAL, MUTA, MUTB) \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const BFIELD& B)                          \
{                                                                                              \
    return Foam::dgCellTmp::binary<OUTFIELD>(A, B, EVAL());                                   \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const BFIELD& B)                    \
{                                                                                              \
    return Foam::dgCellTmp::binary<OUTFIELD>(tA, B, EVAL(), MUTA());                          \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const tmp<BFIELD>& tB)                    \
{                                                                                              \
    return Foam::dgCellTmp::binary<OUTFIELD>(A, tB, EVAL(), MUTB());                          \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const tmp<BFIELD>& tB)              \
{                                                                                              \
    return Foam::dgCellTmp::binary<OUTFIELD>(tA, tB, EVAL(), MUTA(), MUTB());                 \
}

/** \brief Generate unary function overloads for non-tmp/tmp inputs. */
#define DG_CELL_TMP_DEFINE_UNARY_FUNCTION(FUNC, OUTFIELD, INFIELD, EVAL, MUT) \
inline tmp<OUTFIELD> FUNC (const INFIELD& A)                                   \
{                                                                               \
    return Foam::dgCellTmp::unary<OUTFIELD>(A, EVAL());                        \
}                                                                               \
inline tmp<OUTFIELD> FUNC (const tmp<INFIELD>& tA)                             \
{                                                                               \
    return Foam::dgCellTmp::unary<OUTFIELD>(tA, EVAL(), MUT());                \
}

/** \brief Generate `(field OP scalar)` overloads with tmp-aware path. */
#define DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR(OP, OUTFIELD, AFIELD, SCALART, EVAL, MUTA) \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const SCALART& b)                          \
{                                                                                               \
    return Foam::dgCellTmp::binaryFieldScalar<OUTFIELD>(A, b, EVAL());                         \
}                                                                                               \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const SCALART& b)                     \
{                                                                                               \
    return Foam::dgCellTmp::binaryFieldScalar<OUTFIELD>(tA, b, EVAL(), MUTA());                \
}

/** \brief Generate `(scalar OP field)` overloads with tmp-aware path. */
#define DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR(OP, OUTFIELD, SCALART, BFIELD, EVAL, MUTB) \
inline tmp<OUTFIELD> operator OP (const SCALART& a, const BFIELD& B)                         \
{                                                                                              \
    return Foam::dgCellTmp::binaryScalarField<OUTFIELD>(a, B, EVAL());                        \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const SCALART& a, const tmp<BFIELD>& tB)                    \
{                                                                                              \
    return Foam::dgCellTmp::binaryScalarField<OUTFIELD>(a, tB, EVAL(), MUTB());               \
}

/** \brief Generate `(FUNC(field, scalar))` overloads with tmp-aware path. */
#define DG_CELL_TMP_DEFINE_FIELD_SCALAR_FUNCTION(FUNC, OUTFIELD, AFIELD, SCALART, EVAL, MUTA) \
inline tmp<OUTFIELD> FUNC (const AFIELD& A, const SCALART& b)                                   \
{                                                                                                 \
    return Foam::dgCellTmp::binaryFieldScalar<OUTFIELD>(A, b, EVAL());                           \
}                                                                                                 \
inline tmp<OUTFIELD> FUNC (const tmp<AFIELD>& tA, const SCALART& b)                              \
{                                                                                                 \
    return Foam::dgCellTmp::binaryFieldScalar<OUTFIELD>(tA, b, EVAL(), MUTA());                  \
}

/** \brief Generate `(field OP field)` comparison overloads with tmp-aware path. */
#define DG_CELL_TMP_DEFINE_COMPARE_OPERATOR(OP, AFIELD, BFIELD, PRED) \
inline bool operator OP (const AFIELD& A, const BFIELD& B)            \
{                                                                      \
    return Foam::dgCellTmp::compare(A, B, PRED());                     \
}                                                                      \
inline bool operator OP (const tmp<AFIELD>& tA, const BFIELD& B)      \
{                                                                      \
    return Foam::dgCellTmp::compare(tA, B, PRED());                    \
}                                                                      \
inline bool operator OP (const AFIELD& A, const tmp<BFIELD>& tB)      \
{                                                                      \
    return Foam::dgCellTmp::compare(A, tB, PRED());                    \
}                                                                      \
inline bool operator OP (const tmp<AFIELD>& tA, const tmp<BFIELD>& tB)\
{                                                                      \
    return Foam::dgCellTmp::compare(tA, tB, PRED());                   \
}

/** \brief Generate `(field OP scalar)` comparison overloads with tmp-aware path. */
#define DG_CELL_TMP_DEFINE_COMPARE_FIELD_SCALAR(OP, AFIELD, SCALART, PRED) \
inline bool operator OP (const AFIELD& A, const SCALART& b)                \
{                                                                           \
    return Foam::dgCellTmp::compareFieldScalar(A, b, PRED());               \
}                                                                           \
inline bool operator OP (const tmp<AFIELD>& tA, const SCALART& b)          \
{                                                                           \
    return Foam::dgCellTmp::compareFieldScalar(tA, b, PRED());              \
}

/** \brief Generate `(scalar OP field)` comparison overloads with tmp-aware path. */
#define DG_CELL_TMP_DEFINE_COMPARE_SCALAR_FIELD(OP, SCALART, BFIELD, PRED) \
inline bool operator OP (const SCALART& a, const BFIELD& B)                \
{                                                                           \
    return Foam::dgCellTmp::compareScalarField(a, B, PRED());               \
}                                                                           \
inline bool operator OP (const SCALART& a, const tmp<BFIELD>& tB)          \
{                                                                           \
    return Foam::dgCellTmp::compareScalarField(a, tB, PRED());              \
}

#endif
