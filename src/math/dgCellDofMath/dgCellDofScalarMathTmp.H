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
    Mathematical operations between cellDof objects (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofScalarMathTmp_H
#define Foam_dgCellDofScalarMathTmp_H

#include "cellDof.H"
#include "dgCellDofTmpOpCore.H"

namespace Foam
{

struct DgCellDofScaEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

struct DgCellDofScaEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

struct DgCellDofScaEvalMul
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a * b)
    {
        return a * b;
    }
};

struct DgCellDofScaEvalDiv
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a / b)
    {
        return a / b;
    }
};

struct DgCellDofScaMutAddA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a += b; }
};

struct DgCellDofScaMutAddB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b += a; }
};

struct DgCellDofScaMutSubA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a -= b; }
};

struct DgCellDofScaMutSubB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

struct DgCellDofScaMutMulA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a *= b; }
};

struct DgCellDofScaMutMulB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b *= a; }
};

struct DgCellDofScaMutDivA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a /= b; }
};

struct DgCellDofScaMutDivB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a / b; }
};

struct DgCellDofScaPredLess
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

struct DgCellDofScaPredLessEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

struct DgCellDofScaPredGreater
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

struct DgCellDofScaPredGreaterEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

struct DgCellDofScaPredEq
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

struct DgCellDofScaPredNe
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

// Scalar - Scalar
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    +,
    cellDof<scalar>,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaEvalAdd,
    DgCellDofScaMutAddA,
    DgCellDofScaMutAddB
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    -,
    cellDof<scalar>,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaEvalSub,
    DgCellDofScaMutSubA,
    DgCellDofScaMutSubB
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<scalar>,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaEvalMul,
    DgCellDofScaMutMulA,
    DgCellDofScaMutMulB
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellDof<scalar>,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaEvalDiv,
    DgCellDofScaMutDivA,
    DgCellDofScaMutDivB
)

// Scalar - Vector
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<vector>,
    cellDof<scalar>,
    cellDof<vector>,
    DgCellDofScaEvalMul,
    Foam::dgCellDofTmp::NoMutateLeft,
    DgCellDofScaMutMulB
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<vector>,
    cellDof<vector>,
    cellDof<scalar>,
    DgCellDofScaEvalMul,
    DgCellDofScaMutMulA,
    Foam::dgCellDofTmp::NoMutateRight
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellDof<vector>,
    cellDof<vector>,
    cellDof<scalar>,
    DgCellDofScaEvalDiv,
    DgCellDofScaMutDivA,
    Foam::dgCellDofTmp::NoMutateRight
)

// Scalar - tensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<tensor>,
    cellDof<scalar>,
    cellDof<tensor>,
    DgCellDofScaEvalMul,
    Foam::dgCellDofTmp::NoMutateLeft,
    DgCellDofScaMutMulB
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<tensor>,
    cellDof<tensor>,
    cellDof<scalar>,
    DgCellDofScaEvalMul,
    DgCellDofScaMutMulA,
    Foam::dgCellDofTmp::NoMutateRight
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellDof<tensor>,
    cellDof<tensor>,
    cellDof<scalar>,
    DgCellDofScaEvalDiv,
    DgCellDofScaMutDivA,
    Foam::dgCellDofTmp::NoMutateRight
)

// Scalar - symmTensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<symmTensor>,
    cellDof<scalar>,
    cellDof<symmTensor>,
    DgCellDofScaEvalMul,
    Foam::dgCellDofTmp::NoMutateLeft,
    DgCellDofScaMutMulB
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<symmTensor>,
    cellDof<symmTensor>,
    cellDof<scalar>,
    DgCellDofScaEvalMul,
    DgCellDofScaMutMulA,
    Foam::dgCellDofTmp::NoMutateRight
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellDof<symmTensor>,
    cellDof<symmTensor>,
    cellDof<scalar>,
    DgCellDofScaEvalDiv,
    DgCellDofScaMutDivA,
    Foam::dgCellDofTmp::NoMutateRight
)

// Scalar - sphericalTensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<sphericalTensor>,
    cellDof<scalar>,
    cellDof<sphericalTensor>,
    DgCellDofScaEvalMul,
    Foam::dgCellDofTmp::NoMutateLeft,
    DgCellDofScaMutMulB
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    *,
    cellDof<sphericalTensor>,
    cellDof<sphericalTensor>,
    cellDof<scalar>,
    DgCellDofScaEvalMul,
    DgCellDofScaMutMulA,
    Foam::dgCellDofTmp::NoMutateRight
)

DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR
(
    /,
    cellDof<sphericalTensor>,
    cellDof<sphericalTensor>,
    cellDof<scalar>,
    DgCellDofScaEvalDiv,
    DgCellDofScaMutDivA,
    Foam::dgCellDofTmp::NoMutateRight
)

// Comparison operators between cellDof<scalar>
DG_CELL_DOF_TMP_DEFINE_COMPARE_OPERATOR
(
    <,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaPredLess
)

DG_CELL_DOF_TMP_DEFINE_COMPARE_OPERATOR
(
    <=,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaPredLessEq
)

DG_CELL_DOF_TMP_DEFINE_COMPARE_OPERATOR
(
    >,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaPredGreater
)

DG_CELL_DOF_TMP_DEFINE_COMPARE_OPERATOR
(
    >=,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaPredGreaterEq
)

DG_CELL_DOF_TMP_DEFINE_COMPARE_OPERATOR
(
    ==,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaPredEq
)

DG_CELL_DOF_TMP_DEFINE_COMPARE_OPERATOR
(
    !=,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofScaPredNe
)

} // End namespace Foam

#endif

// ************************************************************************* //

