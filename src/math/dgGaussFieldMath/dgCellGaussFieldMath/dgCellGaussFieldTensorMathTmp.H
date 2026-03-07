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
    Tensor operations between cellGaussField objects (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldTensorMathTmp_H
#define Foam_dgCellGaussFieldTensorMathTmp_H

#include "cellGaussField.H"
#include "Tensor.H"
#include "dgCellGaussFieldTmpOpCore.H"

namespace Foam
{

/** \brief Evaluate `dev(.)` at one Gauss point. */
struct DgCellEvalDev
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(dev(a)) { return dev(a); }
};

/** \brief Evaluate `symm(tensor)` at one Gauss point. */
struct DgCellEvalSymm
{
    inline symmTensor operator()(const tensor& a) const { return symm(a); }
};

/** \brief Evaluate `skew(tensor)` at one Gauss point. */
struct DgCellEvalSkew
{
    inline tensor operator()(const tensor& a) const { return skew(a); }
};

/** \brief Evaluate `tr(.)` at one Gauss point. */
struct DgCellEvalTr
{
    template<class A>
    inline scalar operator()(const A& a) const { return tr(a); }
};

/** \brief Evaluate `det(tensor)` at one Gauss point. */
struct DgCellEvalDet
{
    inline scalar operator()(const tensor& a) const { return det(a); }
};

/** \brief Evaluate transpose `T(tensor)` at one Gauss point. */
struct DgCellEvalT
{
    inline tensor operator()(const tensor& a) const { return T(a); }
};

/** \brief Identity mapping for symmTensor transpose overload. */
struct DgCellEvalIdentitySymm
{
    inline symmTensor operator()(const symmTensor& a) const { return a; }
};

/** \brief Evaluate `inv(.)` at one Gauss point. */
struct DgCellEvalInv
{
    template<class A>
    inline tensor operator()(const A& a) const { return inv(a); }
};

/** \brief In-place mutate tensor by `dev(.)`. */
struct DgCellMutDevTensor
{
    inline void operator()(tensor& a) const { a = dev(a); }
};

/** \brief In-place mutate tensor by `skew(.)`. */
struct DgCellMutSkewTensor
{
    inline void operator()(tensor& a) const { a = skew(a); }
};

/** \brief In-place mutate tensor by transpose `T(.)`. */
struct DgCellMutTTensor
{
    inline void operator()(tensor& a) const { a = T(a); }
};

/** \brief In-place mutate symmTensor by `dev(.)`. */
struct DgCellMutDevSymm
{
    inline void operator()(symmTensor& a) const { a = dev(a); }
};

/** \brief No-op mutator for non-reusable tmp paths. */
struct DgCellMutNoopSymm
{
    inline void operator()(symmTensor&) const {}
};

/** \brief Evaluate binary addition at one Gauss point. */
struct DgCellEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

/** \brief Evaluate binary subtraction at one Gauss point. */
struct DgCellEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

/** \brief Evaluate single contraction `&` at one Gauss point. */
struct DgCellEvalAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b)
    {
        return a & b;
    }
};

/** \brief Evaluate double contraction `&&` at one Gauss point. */
struct DgCellEvalAndAnd
{
    template<class A, class B>
    inline scalar operator()(const A& a, const B& b) const
    {
        return a && b;
    }
};

/** \brief In-place mutate left operand for addition. */
struct DgCellMutAddA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a + b; }
};

/** \brief In-place mutate right operand for addition. */
struct DgCellMutAddB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a + b; }
};

/** \brief In-place mutate left operand for subtraction. */
struct DgCellMutSubA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a - b; }
};

/** \brief In-place mutate right operand for subtraction. */
struct DgCellMutSubB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

/** \brief In-place mutate left operand for contraction `&`. */
struct DgCellMutAndA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a & b; }
};

/** \brief In-place mutate right operand for contraction `&`. */
struct DgCellMutAndB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a & b; }
};

/** \brief Unary tensor utilities: `dev`, `symm`, `skew`, `tr`, `det`, `T`, `inv`. */
// Tensor utilities
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(dev,  cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalDev,       DgCellMutDevTensor)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(symm, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalSymm,      Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(skew, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalSkew,      DgCellMutSkewTensor)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(tr,   cellGaussField<scalar>,     cellGaussField<tensor>,     DgCellEvalTr,        Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(det,  cellGaussField<scalar>,     cellGaussField<tensor>,     DgCellEvalDet,       Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(T,    cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalT,         DgCellMutTTensor)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(inv,  cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalInv,       Foam::dgCellTmp::NoMutateUnary)

/** \brief Unary symmTensor utilities: `tr`, `dev`, `T`, `inv`. */
// symmTensor utilities
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(tr,   cellGaussField<scalar>,     cellGaussField<symmTensor>, DgCellEvalTr,        Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(dev,  cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalDev,       DgCellMutDevSymm)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(T,    cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalIdentitySymm, DgCellMutNoopSymm)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(inv,  cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalInv,       Foam::dgCellTmp::NoMutateUnary)

/** \brief Binary tensor-tensor operators: `+`, `-`, `&`, `&&`. */
// tensor - tensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalAdd,    DgCellMutAddA, DgCellMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalSub,    DgCellMutSubA, DgCellMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalAnd,    DgCellMutAndA, DgCellMutAndB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

/** \brief Binary tensor-symmTensor operators: `+`, `-`, `&`, `&&`. */
// tensor - symmTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalAdd,    DgCellMutAddA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalSub,    DgCellMutSubA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalAnd,    DgCellMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

/** \brief Binary symmTensor-tensor operators: `+`, `-`, `&`, `&&`. */
// symmTensor - tensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<tensor>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalAdd,    Foam::dgCellTmp::NoMutateLeft, DgCellMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<tensor>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalSub,    Foam::dgCellTmp::NoMutateLeft, DgCellMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalAnd,    Foam::dgCellTmp::NoMutateLeft, DgCellMutAndB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

/** \brief Binary symmTensor-symmTensor operators: `+`, `-`, `&`, `&&`. */
// symmTensor - symmTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<symmTensor>, cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalAdd,    DgCellMutAddA, DgCellMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<symmTensor>, cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalSub,    DgCellMutSubA, DgCellMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>,     cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalAnd,    Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>,     cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

} // End namespace Foam

#endif

// ************************************************************************* //
