/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Mathematical operations between cellGaussField objects (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldScalarMathTmp_H
#define Foam_dgCellGaussFieldScalarMathTmp_H

#include "cellGaussField.H"
#include "dgCellGaussFieldTmpOpCore.H"

namespace Foam
{

struct DgEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

struct DgEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

struct DgEvalMul
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b)
    {
        return a * b;
    }
};

struct DgEvalDiv
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a / b)
    {
        return a / b;
    }
};

struct DgMutAddA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a += b; }
};

struct DgMutAddB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b += a; }
};

struct DgMutSubA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a -= b; }
};

struct DgMutSubB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

struct DgMutMulA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a *= b; }
};

struct DgMutMulB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b *= a; }
};

struct DgMutDivA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a /= b; }
};

struct DgMutDivB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a / b; }
};

struct DgPredLess
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

struct DgPredLessEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

struct DgPredGreater
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

struct DgPredGreaterEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

struct DgPredEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

struct DgPredNe
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

// Scalar - Scalar
DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    +,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgEvalAdd,
    DgMutAddA,
    DgMutAddB
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    -,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgEvalSub,
    DgMutSubA,
    DgMutSubB
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgEvalMul,
    DgMutMulA,
    DgMutMulB
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgEvalDiv,
    DgMutDivA,
    DgMutDivB
)

// Scalar - Vector
DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<vector>,
    cellGaussField<scalar>,
    cellGaussField<vector>,
    DgEvalMul,
    Foam::dgCellTmp::NoMutateLeft,
    DgMutMulB
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<vector>,
    cellGaussField<vector>,
    cellGaussField<scalar>,
    DgEvalMul,
    DgMutMulA,
    Foam::dgCellTmp::NoMutateRight
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellGaussField<vector>,
    cellGaussField<vector>,
    cellGaussField<scalar>,
    DgEvalDiv,
    DgMutDivA,
    Foam::dgCellTmp::NoMutateRight
)

// Scalar - tensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<tensor>,
    cellGaussField<scalar>,
    cellGaussField<tensor>,
    DgEvalMul,
    Foam::dgCellTmp::NoMutateLeft,
    DgMutMulB
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<tensor>,
    cellGaussField<tensor>,
    cellGaussField<scalar>,
    DgEvalMul,
    DgMutMulA,
    Foam::dgCellTmp::NoMutateRight
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellGaussField<tensor>,
    cellGaussField<tensor>,
    cellGaussField<scalar>,
    DgEvalDiv,
    DgMutDivA,
    Foam::dgCellTmp::NoMutateRight
)

// Scalar - symmTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<symmTensor>,
    cellGaussField<scalar>,
    cellGaussField<symmTensor>,
    DgEvalMul,
    Foam::dgCellTmp::NoMutateLeft,
    DgMutMulB
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<symmTensor>,
    cellGaussField<symmTensor>,
    cellGaussField<scalar>,
    DgEvalMul,
    DgMutMulA,
    Foam::dgCellTmp::NoMutateRight
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellGaussField<symmTensor>,
    cellGaussField<symmTensor>,
    cellGaussField<scalar>,
    DgEvalDiv,
    DgMutDivA,
    Foam::dgCellTmp::NoMutateRight
)

// Scalar - sphericalTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<sphericalTensor>,
    cellGaussField<scalar>,
    cellGaussField<sphericalTensor>,
    DgEvalMul,
    Foam::dgCellTmp::NoMutateLeft,
    DgMutMulB
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellGaussField<sphericalTensor>,
    cellGaussField<sphericalTensor>,
    cellGaussField<scalar>,
    DgEvalMul,
    DgMutMulA,
    Foam::dgCellTmp::NoMutateRight
)

DG_CELL_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellGaussField<sphericalTensor>,
    cellGaussField<sphericalTensor>,
    cellGaussField<scalar>,
    DgEvalDiv,
    DgMutDivA,
    Foam::dgCellTmp::NoMutateRight
)

// Comparison operators between cellGaussField<scalar>
DG_CELL_TMP_DEFINE_COMPARE_OPERATOR
(
    <,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgPredLess
)

DG_CELL_TMP_DEFINE_COMPARE_OPERATOR
(
    <=,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgPredLessEq
)

DG_CELL_TMP_DEFINE_COMPARE_OPERATOR
(
    >,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgPredGreater
)

DG_CELL_TMP_DEFINE_COMPARE_OPERATOR
(
    >=,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgPredGreaterEq
)

DG_CELL_TMP_DEFINE_COMPARE_OPERATOR
(
    ==,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgPredEq
)

DG_CELL_TMP_DEFINE_COMPARE_OPERATOR
(
    !=,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgPredNe
)

} // End namespace Foam

#endif

// ************************************************************************* //

