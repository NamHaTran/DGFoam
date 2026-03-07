/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Core helpers/macros to reduce tmp/non-tmp overload boilerplate for
    faceGaussField-style plus/minus Gauss-point storage.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldTmpOpCore_H
#define Foam_dgFaceGaussFieldTmpOpCore_H

#include "tmp.H"
#include <type_traits>

namespace Foam
{
namespace dgFaceTmp
{

struct NoMutateLeft
{
    template<class LhsElem, class RhsElem>
    inline void operator()
    (
        LhsElem&,
        LhsElem&,
        const RhsElem&,
        const RhsElem&
    ) const {}
};

struct NoMutateRight
{
    template<class LhsElem, class RhsElem>
    inline void operator()
    (
        const LhsElem&,
        const LhsElem&,
        RhsElem&,
        RhsElem&
    ) const {}
};

struct NoMutateUnary
{
    template<class Elem>
    inline void operator()(Elem&, Elem&) const {}
};

template<class OutField, class RefField>
inline tmp<OutField> newLike(const RefField& A)
{
    return tmp<OutField>(new OutField(A.cellID(), A.dgMesh()));
}

template<class OutField, class InField, class Eval>
inline tmp<OutField> unary
(
    const InField& A,
    const Eval& eval
)
{
    auto tRes = newLike<OutField>(A);
    OutField& R = tRes.ref();

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        R.plusValueAt(gp) = eval(A.plusValue(gp));
        R.minusValueAt(gp) = eval(A.minusValue(gp));
    }

    return tRes;
}

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

        for (label gp = 0; gp < A.nGauss(); ++gp)
        {
            mutate(A.plusValueAt(gp), A.minusValueAt(gp));
        }

        return tA;
    }

    return unary<OutField>(tA(), eval);
}

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

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        R.plusValueAt(gp) = eval(A.plusValue(gp), B.plusValue(gp));
        R.minusValueAt(gp) = eval(A.minusValue(gp), B.minusValue(gp));
    }

    return tRes;
}

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

        for (label gp = 0; gp < A.nGauss(); ++gp)
        {
            mutateA
            (
                A.plusValueAt(gp),
                A.minusValueAt(gp),
                B.plusValue(gp),
                B.minusValue(gp)
            );
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

        for (label gp = 0; gp < B.nGauss(); ++gp)
        {
            mutateB
            (
                A.plusValue(gp),
                A.minusValue(gp),
                B.plusValueAt(gp),
                B.minusValueAt(gp)
            );
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

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        R.plusValueAt(gp) = eval(A.plusValue(gp), b);
        R.minusValueAt(gp) = eval(A.minusValue(gp), b);
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

        for (label gp = 0; gp < A.nGauss(); ++gp)
        {
            mutateA(A.plusValueAt(gp), A.minusValueAt(gp), b);
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

    for (label gp = 0; gp < B.nGauss(); ++gp)
    {
        R.plusValueAt(gp) = eval(a, B.plusValue(gp));
        R.minusValueAt(gp) = eval(a, B.minusValue(gp));
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

        for (label gp = 0; gp < B.nGauss(); ++gp)
        {
            mutateB(a, B.plusValueAt(gp), B.minusValueAt(gp));
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

template<class AField, class BField, class Pred>
inline bool compare
(
    const AField& A,
    const BField& B,
    const Pred& pred
)
{
    if (A.nGauss() != B.nGauss())
    {
        return false;
    }

    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        if (!pred(A.plusValue(gp), B.plusValue(gp)))
        {
            return false;
        }

        if (!pred(A.minusValue(gp), B.minusValue(gp)))
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

template<class AField, class ScalarType, class Pred>
inline bool compareFieldScalar
(
    const AField& A,
    const ScalarType& b,
    const Pred& pred
)
{
    for (label gp = 0; gp < A.nGauss(); ++gp)
    {
        if (!pred(A.plusValue(gp), b))
        {
            return false;
        }

        if (!pred(A.minusValue(gp), b))
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

template<class ScalarType, class BField, class Pred>
inline bool compareScalarField
(
    const ScalarType& a,
    const BField& B,
    const Pred& pred
)
{
    for (label gp = 0; gp < B.nGauss(); ++gp)
    {
        if (!pred(a, B.plusValue(gp)))
        {
            return false;
        }

        if (!pred(a, B.minusValue(gp)))
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

} // End namespace dgFaceTmp
} // End namespace Foam

#define DG_FACE_TMP_DEFINE_BINARY_OPERATOR(OP, OUTFIELD, AFIELD, BFIELD, EVAL, MUTA, MUTB) \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const BFIELD& B)                          \
{                                                                                              \
    return Foam::dgFaceTmp::binary<OUTFIELD>(A, B, EVAL());                                   \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const BFIELD& B)                    \
{                                                                                              \
    return Foam::dgFaceTmp::binary<OUTFIELD>(tA, B, EVAL(), MUTA());                          \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const tmp<BFIELD>& tB)                    \
{                                                                                              \
    return Foam::dgFaceTmp::binary<OUTFIELD>(A, tB, EVAL(), MUTB());                          \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const tmp<BFIELD>& tB)              \
{                                                                                              \
    return Foam::dgFaceTmp::binary<OUTFIELD>(tA, tB, EVAL(), MUTA(), MUTB());                 \
}

#define DG_FACE_TMP_DEFINE_UNARY_FUNCTION(FUNC, OUTFIELD, INFIELD, EVAL, MUT) \
inline tmp<OUTFIELD> FUNC (const INFIELD& A)                                   \
{                                                                               \
    return Foam::dgFaceTmp::unary<OUTFIELD>(A, EVAL());                        \
}                                                                               \
inline tmp<OUTFIELD> FUNC (const tmp<INFIELD>& tA)                             \
{                                                                               \
    return Foam::dgFaceTmp::unary<OUTFIELD>(tA, EVAL(), MUT());                \
}

#define DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR(OP, OUTFIELD, AFIELD, SCALART, EVAL, MUTA) \
inline tmp<OUTFIELD> operator OP (const AFIELD& A, const SCALART& b)                          \
{                                                                                               \
    return Foam::dgFaceTmp::binaryFieldScalar<OUTFIELD>(A, b, EVAL());                         \
}                                                                                               \
inline tmp<OUTFIELD> operator OP (const tmp<AFIELD>& tA, const SCALART& b)                     \
{                                                                                               \
    return Foam::dgFaceTmp::binaryFieldScalar<OUTFIELD>(tA, b, EVAL(), MUTA());                \
}

#define DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR(OP, OUTFIELD, SCALART, BFIELD, EVAL, MUTB) \
inline tmp<OUTFIELD> operator OP (const SCALART& a, const BFIELD& B)                         \
{                                                                                              \
    return Foam::dgFaceTmp::binaryScalarField<OUTFIELD>(a, B, EVAL());                        \
}                                                                                              \
inline tmp<OUTFIELD> operator OP (const SCALART& a, const tmp<BFIELD>& tB)                    \
{                                                                                              \
    return Foam::dgFaceTmp::binaryScalarField<OUTFIELD>(a, tB, EVAL(), MUTB());               \
}

#define DG_FACE_TMP_DEFINE_FIELD_SCALAR_FUNCTION(FUNC, OUTFIELD, AFIELD, SCALART, EVAL, MUTA) \
inline tmp<OUTFIELD> FUNC (const AFIELD& A, const SCALART& b)                                   \
{                                                                                                 \
    return Foam::dgFaceTmp::binaryFieldScalar<OUTFIELD>(A, b, EVAL());                           \
}                                                                                                 \
inline tmp<OUTFIELD> FUNC (const tmp<AFIELD>& tA, const SCALART& b)                              \
{                                                                                                 \
    return Foam::dgFaceTmp::binaryFieldScalar<OUTFIELD>(tA, b, EVAL(), MUTA());                  \
}

#define DG_FACE_TMP_DEFINE_COMPARE_OPERATOR(OP, AFIELD, BFIELD, PRED) \
inline bool operator OP (const AFIELD& A, const BFIELD& B)            \
{                                                                      \
    return Foam::dgFaceTmp::compare(A, B, PRED());                     \
}                                                                      \
inline bool operator OP (const tmp<AFIELD>& tA, const BFIELD& B)      \
{                                                                      \
    return Foam::dgFaceTmp::compare(tA, B, PRED());                    \
}                                                                      \
inline bool operator OP (const AFIELD& A, const tmp<BFIELD>& tB)      \
{                                                                      \
    return Foam::dgFaceTmp::compare(A, tB, PRED());                    \
}                                                                      \
inline bool operator OP (const tmp<AFIELD>& tA, const tmp<BFIELD>& tB)\
{                                                                      \
    return Foam::dgFaceTmp::compare(tA, tB, PRED());                   \
}

#define DG_FACE_TMP_DEFINE_COMPARE_FIELD_SCALAR(OP, AFIELD, SCALART, PRED) \
inline bool operator OP (const AFIELD& A, const SCALART& b)                \
{                                                                           \
    return Foam::dgFaceTmp::compareFieldScalar(A, b, PRED());               \
}                                                                           \
inline bool operator OP (const tmp<AFIELD>& tA, const SCALART& b)          \
{                                                                           \
    return Foam::dgFaceTmp::compareFieldScalar(tA, b, PRED());              \
}

#define DG_FACE_TMP_DEFINE_COMPARE_SCALAR_FIELD(OP, SCALART, BFIELD, PRED) \
inline bool operator OP (const SCALART& a, const BFIELD& B)                \
{                                                                           \
    return Foam::dgFaceTmp::compareScalarField(a, B, PRED());               \
}                                                                           \
inline bool operator OP (const SCALART& a, const tmp<BFIELD>& tB)          \
{                                                                           \
    return Foam::dgFaceTmp::compareScalarField(a, tB, PRED());              \
}

#endif
