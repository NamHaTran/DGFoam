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
    Mathematical operations between cellGaussField objects (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldScalarMathTmp_H
#define Foam_dgCellGaussFieldScalarMathTmp_H

#include "cellGaussField.H"
#include "dgCellGaussFieldTmpOpCore.H"

namespace Foam
{

/** \brief Evaluate binary addition at one Gauss point. */
struct DgEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

/** \brief Evaluate binary subtraction at one Gauss point. */
struct DgEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

/** \brief Evaluate binary multiplication at one Gauss point. */
struct DgEvalMul
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b)
    {
        return a * b;
    }
};

/** \brief Evaluate binary division at one Gauss point. */
struct DgEvalDiv
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a / b)
    {
        return a / b;
    }
};

/** \brief In-place mutate left operand for addition. */
struct DgMutAddA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a += b; }
};

/** \brief In-place mutate right operand for addition. */
struct DgMutAddB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b += a; }
};

/** \brief In-place mutate left operand for subtraction. */
struct DgMutSubA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a -= b; }
};

/** \brief In-place mutate right operand for subtraction. */
struct DgMutSubB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

/** \brief In-place mutate left operand for multiplication. */
struct DgMutMulA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a *= b; }
};

/** \brief In-place mutate right operand for multiplication. */
struct DgMutMulB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b *= a; }
};

/** \brief In-place mutate left operand for division. */
struct DgMutDivA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a /= b; }
};

/** \brief In-place mutate right operand for division. */
struct DgMutDivB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a / b; }
};

/** \brief Predicate functor for `<` on scalar fields. */
struct DgPredLess
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

/** \brief Predicate functor for `<=` on scalar fields. */
struct DgPredLessEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

/** \brief Predicate functor for `>` on scalar fields. */
struct DgPredGreater
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

/** \brief Predicate functor for `>=` on scalar fields. */
struct DgPredGreaterEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

/** \brief Predicate functor for `==` on scalar fields. */
struct DgPredEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

/** \brief Predicate functor for `!=` on scalar fields. */
struct DgPredNe
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

/** \brief Scalar-scalar operators: `+`, `-`, `*`, `/`. */
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

/** \brief Scalar-vector operators: `*` and vector `/` scalar. */
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

/** \brief Scalar-tensor operators: `*` and tensor `/` scalar. */
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

/** \brief Scalar-symmTensor operators: `*` and symmTensor `/` scalar. */
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

/** \brief Scalar-sphericalTensor operators: `*` and sphericalTensor `/` scalar. */
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

/** \brief Comparison operators between two scalar cellGaussField objects. */
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
