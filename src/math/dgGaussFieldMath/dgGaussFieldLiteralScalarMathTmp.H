/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Operator overloads for arithmetic between literal scalar and GaussField
    with tmp-aware implementations.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldLiteralScalarMathTmp_H
#define Foam_dgGaussFieldLiteralScalarMathTmp_H

#include "GaussField.H"
#include "dgFaceGaussFieldLiteralScalarMathTmp.H"
#include "dgCellGaussFieldLiteralScalarMathTmp.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"

namespace Foam
{

/** \brief Allocate output GaussField with geometry copied from reference field. */
template<class Type>
inline tmp<GaussField<Type>> dgGaussNewLike(const GaussField<Type>& A)
{
    return tmp<GaussField<Type>>(new GaussField<Type>(A.cellID(), A.dgMesh()));
}

/** \brief Multiply GaussField by scalar literal (`s * A`). */
template<class Type>
inline tmp<GaussField<Type>> dgMul
(
    const scalar& s,
    const GaussField<Type>& A
)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<Type>& R = tRes.ref();
    R.cellField() = s * A.cellField();
    R.faceField() = s * A.faceField();
    return tRes;
}

/** \brief Multiply tmp GaussField by scalar literal with in-place reuse when possible. */
template<class Type>
inline tmp<GaussField<Type>> dgMul
(
    const scalar& s,
    const tmp<GaussField<Type>>& tA
)
{
    if (tA.isTmp())
    {
        GaussField<Type>& A = const_cast<tmp<GaussField<Type>>&>(tA).ref();
        A.cellField() = s * A.cellField();
        A.faceField() = s * A.faceField();
        return tA;
    }

    return dgMul(s, tA());
}

/** \brief Divide GaussField by scalar literal (`A / s`). */
template<class Type>
inline tmp<GaussField<Type>> dgDiv
(
    const GaussField<Type>& A,
    const scalar& s
)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<Type>& R = tRes.ref();
    R.cellField() = A.cellField() / s;
    R.faceField() = A.faceField() / s;
    return tRes;
}

/** \brief Divide tmp GaussField by scalar literal with in-place reuse when possible. */
template<class Type>
inline tmp<GaussField<Type>> dgDiv
(
    const tmp<GaussField<Type>>& tA,
    const scalar& s
)
{
    if (tA.isTmp())
    {
        GaussField<Type>& A = const_cast<tmp<GaussField<Type>>&>(tA).ref();
        A.cellField() = A.cellField() / s;
        A.faceField() = A.faceField() / s;
        return tA;
    }

    return dgDiv(tA(), s);
}

/** \brief Left division for scalar field (`s / A`). */
inline tmp<GaussField<scalar>> dgDivLeft
(
    const scalar& s,
    const GaussField<scalar>& A
)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<scalar>& R = tRes.ref();
    R.cellField() = s / A.cellField();
    R.faceField() = s / A.faceField();
    return tRes;
}

/** \brief Left division for tmp scalar field (`s / tA`). */
inline tmp<GaussField<scalar>> dgDivLeft
(
    const scalar& s,
    const tmp<GaussField<scalar>>& tA
)
{
    if (tA.isTmp())
    {
        GaussField<scalar>& A = const_cast<tmp<GaussField<scalar>>&>(tA).ref();
        A.cellField() = s / A.cellField();
        A.faceField() = s / A.faceField();
        return tA;
    }

    return dgDivLeft(s, tA());
}

/** \brief Add scalar literal to scalar field (`s + A`). */
inline tmp<GaussField<scalar>> dgAdd
(
    const scalar& s,
    const GaussField<scalar>& A
)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<scalar>& R = tRes.ref();
    R.cellField() = s + A.cellField();
    R.faceField() = s + A.faceField();
    return tRes;
}

/** \brief Add scalar literal to tmp scalar field (`s + tA`). */
inline tmp<GaussField<scalar>> dgAdd
(
    const scalar& s,
    const tmp<GaussField<scalar>>& tA
)
{
    if (tA.isTmp())
    {
        GaussField<scalar>& A = const_cast<tmp<GaussField<scalar>>&>(tA).ref();
        A.cellField() = s + A.cellField();
        A.faceField() = s + A.faceField();
        return tA;
    }

    return dgAdd(s, tA());
}

/** \brief Left subtraction for scalar field (`s - A`). */
inline tmp<GaussField<scalar>> dgSubLeft
(
    const scalar& s,
    const GaussField<scalar>& A
)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<scalar>& R = tRes.ref();
    R.cellField() = s - A.cellField();
    R.faceField() = s - A.faceField();
    return tRes;
}

/** \brief Left subtraction for tmp scalar field (`s - tA`). */
inline tmp<GaussField<scalar>> dgSubLeft
(
    const scalar& s,
    const tmp<GaussField<scalar>>& tA
)
{
    if (tA.isTmp())
    {
        GaussField<scalar>& A = const_cast<tmp<GaussField<scalar>>&>(tA).ref();
        A.cellField() = s - A.cellField();
        A.faceField() = s - A.faceField();
        return tA;
    }

    return dgSubLeft(s, tA());
}

/** \brief Right subtraction for scalar field (`A - s`). */
inline tmp<GaussField<scalar>> dgSubRight
(
    const GaussField<scalar>& A,
    const scalar& s
)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<scalar>& R = tRes.ref();
    R.cellField() = A.cellField() - s;
    R.faceField() = A.faceField() - s;
    return tRes;
}

/** \brief Right subtraction for tmp scalar field (`tA - s`). */
inline tmp<GaussField<scalar>> dgSubRight
(
    const tmp<GaussField<scalar>>& tA,
    const scalar& s
)
{
    if (tA.isTmp())
    {
        GaussField<scalar>& A = const_cast<tmp<GaussField<scalar>>&>(tA).ref();
        A.cellField() = A.cellField() - s;
        A.faceField() = A.faceField() - s;
        return tA;
    }

    return dgSubRight(tA(), s);
}

/** \brief Compare all scalar field values against a scalar literal. */
template<class Pred>
inline bool dgCompareFieldScalar
(
    const GaussField<scalar>& A,
    const scalar& b,
    const Pred& pred
)
{
    return pred(A.cellField(), b) && pred(A.faceField(), b);
}

/** \brief Compare tmp scalar field values against a scalar literal. */
template<class Pred>
inline bool dgCompareFieldScalar
(
    const tmp<GaussField<scalar>>& tA,
    const scalar& b,
    const Pred& pred
)
{
    return dgCompareFieldScalar(tA(), b, pred);
}

/** \brief Compare a scalar literal against all scalar field values. */
template<class Pred>
inline bool dgCompareScalarField
(
    const scalar& a,
    const GaussField<scalar>& B,
    const Pred& pred
)
{
    return pred(a, B.cellField()) && pred(a, B.faceField());
}

/** \brief Compare a scalar literal against tmp scalar field values. */
template<class Pred>
inline bool dgCompareScalarField
(
    const scalar& a,
    const tmp<GaussField<scalar>>& tB,
    const Pred& pred
)
{
    return dgCompareScalarField(a, tB(), pred);
}

/** \brief Predicate functor for `<` comparison. */
struct DgGaussLessPred
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a < b; }
};

/** \brief Predicate functor for `<=` comparison. */
struct DgGaussLessEqPred
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a <= b; }
};

/** \brief Predicate functor for `>` comparison. */
struct DgGaussGreaterPred
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a > b; }
};

/** \brief Predicate functor for `>=` comparison. */
struct DgGaussGreaterEqPred
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a >= b; }
};

/** \brief Predicate functor for `==` comparison. */
struct DgGaussEqPred
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a == b; }
};

/** \brief Predicate functor for `!=` comparison. */
struct DgGaussNePred
{
    template<class A, class B>
    inline bool operator()(const A& a, const B& b) const { return a != b; }
};

inline tmp<GaussField<scalar>> operator*(const scalar& s, const GaussField<scalar>& A) { return dgMul(s, A); }
inline tmp<GaussField<scalar>> operator*(const scalar& s, const tmp<GaussField<scalar>>& tA) { return dgMul(s, tA); }
inline tmp<GaussField<scalar>> operator*(const GaussField<scalar>& A, const scalar& s) { return dgMul(s, A); }
inline tmp<GaussField<scalar>> operator*(const tmp<GaussField<scalar>>& tA, const scalar& s) { return dgMul(s, tA); }

inline tmp<GaussField<vector>> operator*(const scalar& s, const GaussField<vector>& A) { return dgMul(s, A); }
inline tmp<GaussField<vector>> operator*(const scalar& s, const tmp<GaussField<vector>>& tA) { return dgMul(s, tA); }
inline tmp<GaussField<vector>> operator*(const GaussField<vector>& A, const scalar& s) { return dgMul(s, A); }
inline tmp<GaussField<vector>> operator*(const tmp<GaussField<vector>>& tA, const scalar& s) { return dgMul(s, tA); }

inline tmp<GaussField<tensor>> operator*(const scalar& s, const GaussField<tensor>& A) { return dgMul(s, A); }
inline tmp<GaussField<tensor>> operator*(const scalar& s, const tmp<GaussField<tensor>>& tA) { return dgMul(s, tA); }
inline tmp<GaussField<tensor>> operator*(const GaussField<tensor>& A, const scalar& s) { return dgMul(s, A); }
inline tmp<GaussField<tensor>> operator*(const tmp<GaussField<tensor>>& tA, const scalar& s) { return dgMul(s, tA); }

inline tmp<GaussField<symmTensor>> operator*(const scalar& s, const GaussField<symmTensor>& A) { return dgMul(s, A); }
inline tmp<GaussField<symmTensor>> operator*(const scalar& s, const tmp<GaussField<symmTensor>>& tA) { return dgMul(s, tA); }
inline tmp<GaussField<symmTensor>> operator*(const GaussField<symmTensor>& A, const scalar& s) { return dgMul(s, A); }
inline tmp<GaussField<symmTensor>> operator*(const tmp<GaussField<symmTensor>>& tA, const scalar& s) { return dgMul(s, tA); }

inline tmp<GaussField<sphericalTensor>> operator*(const scalar& s, const GaussField<sphericalTensor>& A) { return dgMul(s, A); }
inline tmp<GaussField<sphericalTensor>> operator*(const scalar& s, const tmp<GaussField<sphericalTensor>>& tA) { return dgMul(s, tA); }
inline tmp<GaussField<sphericalTensor>> operator*(const GaussField<sphericalTensor>& A, const scalar& s) { return dgMul(s, A); }
inline tmp<GaussField<sphericalTensor>> operator*(const tmp<GaussField<sphericalTensor>>& tA, const scalar& s) { return dgMul(s, tA); }

inline tmp<GaussField<scalar>> operator/(const scalar& s, const GaussField<scalar>& A) { return dgDivLeft(s, A); }
inline tmp<GaussField<scalar>> operator/(const scalar& s, const tmp<GaussField<scalar>>& tA) { return dgDivLeft(s, tA); }
inline tmp<GaussField<scalar>> operator/(const GaussField<scalar>& A, const scalar& s) { return dgDiv(A, s); }
inline tmp<GaussField<scalar>> operator/(const tmp<GaussField<scalar>>& tA, const scalar& s) { return dgDiv(tA, s); }

inline tmp<GaussField<vector>> operator/(const GaussField<vector>& A, const scalar& s) { return dgDiv(A, s); }
inline tmp<GaussField<vector>> operator/(const tmp<GaussField<vector>>& tA, const scalar& s) { return dgDiv(tA, s); }
inline tmp<GaussField<tensor>> operator/(const GaussField<tensor>& A, const scalar& s) { return dgDiv(A, s); }
inline tmp<GaussField<tensor>> operator/(const tmp<GaussField<tensor>>& tA, const scalar& s) { return dgDiv(tA, s); }
inline tmp<GaussField<symmTensor>> operator/(const GaussField<symmTensor>& A, const scalar& s) { return dgDiv(A, s); }
inline tmp<GaussField<symmTensor>> operator/(const tmp<GaussField<symmTensor>>& tA, const scalar& s) { return dgDiv(tA, s); }
inline tmp<GaussField<sphericalTensor>> operator/(const GaussField<sphericalTensor>& A, const scalar& s) { return dgDiv(A, s); }
inline tmp<GaussField<sphericalTensor>> operator/(const tmp<GaussField<sphericalTensor>>& tA, const scalar& s) { return dgDiv(tA, s); }

inline tmp<GaussField<scalar>> operator+(const scalar& s, const GaussField<scalar>& A) { return dgAdd(s, A); }
inline tmp<GaussField<scalar>> operator+(const scalar& s, const tmp<GaussField<scalar>>& tA) { return dgAdd(s, tA); }
inline tmp<GaussField<scalar>> operator+(const GaussField<scalar>& A, const scalar& s) { return dgAdd(s, A); }
inline tmp<GaussField<scalar>> operator+(const tmp<GaussField<scalar>>& tA, const scalar& s) { return dgAdd(s, tA); }

inline tmp<GaussField<scalar>> operator-(const scalar& s, const GaussField<scalar>& A) { return dgSubLeft(s, A); }
inline tmp<GaussField<scalar>> operator-(const scalar& s, const tmp<GaussField<scalar>>& tA) { return dgSubLeft(s, tA); }
inline tmp<GaussField<scalar>> operator-(const GaussField<scalar>& A, const scalar& s) { return dgSubRight(A, s); }
inline tmp<GaussField<scalar>> operator-(const tmp<GaussField<scalar>>& tA, const scalar& s) { return dgSubRight(tA, s); }

inline tmp<GaussField<scalar>> pow(const GaussField<scalar>& A, const scalar& s)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<scalar>& R = tRes.ref();
    R.cellField() = pow(A.cellField(), s);
    R.faceField() = pow(A.faceField(), s);
    return tRes;
}

inline tmp<GaussField<scalar>> pow(const tmp<GaussField<scalar>>& tA, const scalar& s)
{
    if (tA.isTmp())
    {
        GaussField<scalar>& A = const_cast<tmp<GaussField<scalar>>&>(tA).ref();
        A.cellField() = pow(A.cellField(), s);
        A.faceField() = pow(A.faceField(), s);
        return tA;
    }
    return pow(tA(), s);
}

inline tmp<GaussField<scalar>> sqrt(const GaussField<scalar>& A)
{
    auto tRes = dgGaussNewLike(A);
    GaussField<scalar>& R = tRes.ref();
    R.cellField() = sqrt(A.cellField());
    R.faceField() = sqrt(A.faceField());
    return tRes;
}

inline tmp<GaussField<scalar>> sqrt(const tmp<GaussField<scalar>>& tA)
{
    if (tA.isTmp())
    {
        GaussField<scalar>& A = const_cast<tmp<GaussField<scalar>>&>(tA).ref();
        A.cellField() = sqrt(A.cellField());
        A.faceField() = sqrt(A.faceField());
        return tA;
    }
    return sqrt(tA());
}

inline bool operator>(const GaussField<scalar>& A, const scalar b) { return dgCompareFieldScalar(A, b, DgGaussGreaterPred()); }
inline bool operator>(const tmp<GaussField<scalar>>& tA, const scalar b) { return dgCompareFieldScalar(tA, b, DgGaussGreaterPred()); }
inline bool operator>(const scalar a, const GaussField<scalar>& B) { return dgCompareScalarField(a, B, DgGaussGreaterPred()); }
inline bool operator>(const scalar a, const tmp<GaussField<scalar>>& tB) { return dgCompareScalarField(a, tB, DgGaussGreaterPred()); }

inline bool operator>=(const GaussField<scalar>& A, const scalar b) { return dgCompareFieldScalar(A, b, DgGaussGreaterEqPred()); }
inline bool operator>=(const tmp<GaussField<scalar>>& tA, const scalar b) { return dgCompareFieldScalar(tA, b, DgGaussGreaterEqPred()); }
inline bool operator>=(const scalar a, const GaussField<scalar>& B) { return dgCompareScalarField(a, B, DgGaussGreaterEqPred()); }
inline bool operator>=(const scalar a, const tmp<GaussField<scalar>>& tB) { return dgCompareScalarField(a, tB, DgGaussGreaterEqPred()); }

inline bool operator<(const GaussField<scalar>& A, const scalar b) { return dgCompareFieldScalar(A, b, DgGaussLessPred()); }
inline bool operator<(const tmp<GaussField<scalar>>& tA, const scalar b) { return dgCompareFieldScalar(tA, b, DgGaussLessPred()); }
inline bool operator<(const scalar a, const GaussField<scalar>& B) { return dgCompareScalarField(a, B, DgGaussLessPred()); }
inline bool operator<(const scalar a, const tmp<GaussField<scalar>>& tB) { return dgCompareScalarField(a, tB, DgGaussLessPred()); }

inline bool operator<=(const GaussField<scalar>& A, const scalar b) { return dgCompareFieldScalar(A, b, DgGaussLessEqPred()); }
inline bool operator<=(const tmp<GaussField<scalar>>& tA, const scalar b) { return dgCompareFieldScalar(tA, b, DgGaussLessEqPred()); }
inline bool operator<=(const scalar a, const GaussField<scalar>& B) { return dgCompareScalarField(a, B, DgGaussLessEqPred()); }
inline bool operator<=(const scalar a, const tmp<GaussField<scalar>>& tB) { return dgCompareScalarField(a, tB, DgGaussLessEqPred()); }

inline bool operator==(const GaussField<scalar>& A, const scalar b) { return dgCompareFieldScalar(A, b, DgGaussEqPred()); }
inline bool operator==(const tmp<GaussField<scalar>>& tA, const scalar b) { return dgCompareFieldScalar(tA, b, DgGaussEqPred()); }
inline bool operator==(const scalar a, const GaussField<scalar>& B) { return dgCompareScalarField(a, B, DgGaussEqPred()); }
inline bool operator==(const scalar a, const tmp<GaussField<scalar>>& tB) { return dgCompareScalarField(a, tB, DgGaussEqPred()); }

inline bool operator!=(const GaussField<scalar>& A, const scalar b) { return dgCompareFieldScalar(A, b, DgGaussNePred()); }
inline bool operator!=(const tmp<GaussField<scalar>>& tA, const scalar b) { return dgCompareFieldScalar(tA, b, DgGaussNePred()); }
inline bool operator!=(const scalar a, const GaussField<scalar>& B) { return dgCompareScalarField(a, B, DgGaussNePred()); }
inline bool operator!=(const scalar a, const tmp<GaussField<scalar>>& tB) { return dgCompareScalarField(a, tB, DgGaussNePred()); }

} // End namespace Foam

#endif

// ************************************************************************* //
