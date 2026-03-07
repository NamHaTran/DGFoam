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
    Vector math operations between cellGaussField objects (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldVectorMathTmp_H
#define Foam_dgCellGaussFieldVectorMathTmp_H

#include "cellGaussField.H"
#include "VectorSpace.H"
#include "dgCellGaussFieldTmpOpCore.H"

namespace Foam
{

/** \brief Evaluate vector addition at one Gauss point. */
struct DgCellVecEvalAdd
{
    inline vector operator()(const vector& a, const vector& b) const { return a + b; }
};

/** \brief Evaluate vector subtraction at one Gauss point. */
struct DgCellVecEvalSub
{
    inline vector operator()(const vector& a, const vector& b) const { return a - b; }
};

/** \brief Evaluate vector cross product at one Gauss point. */
struct DgCellVecEvalCross
{
    inline vector operator()(const vector& a, const vector& b) const { return a ^ b; }
};

/** \brief Evaluate vector dot product at one Gauss point. */
struct DgCellVecEvalDot
{
    inline scalar operator()(const vector& a, const vector& b) const { return a & b; }
};

/** \brief Evaluate vector outer product at one Gauss point. */
struct DgCellVecEvalOuter
{
    inline tensor operator()(const vector& a, const vector& b) const { return a * b; }
};

/** \brief Evaluate vector-tensor contraction `vector & tensorLike`. */
struct DgCellVecEvalAndVT
{
    template<class T>
    inline vector operator()(const vector& a, const T& b) const { return a & b; }
};

/** \brief Evaluate tensor-vector contraction `tensorLike & vector`. */
struct DgCellVecEvalAndTV
{
    template<class T>
    inline vector operator()(const T& a, const vector& b) const { return a & b; }
};

/** \brief Evaluate vector magnitude at one Gauss point. */
struct DgCellVecEvalMag
{
    inline scalar operator()(const vector& a) const { return mag(a); }
};

/** \brief Evaluate squared vector magnitude at one Gauss point. */
struct DgCellVecEvalMagSqr
{
    inline scalar operator()(const vector& a) const { return magSqr(a); }
};

/** \brief In-place mutate left operand for vector addition. */
struct DgCellVecMutAddA
{
    inline void operator()(vector& a, const vector& b) const { a += b; }
};

/** \brief In-place mutate right operand for vector addition. */
struct DgCellVecMutAddB
{
    inline void operator()(const vector& a, vector& b) const { b += a; }
};

/** \brief In-place mutate left operand for vector subtraction. */
struct DgCellVecMutSubA
{
    inline void operator()(vector& a, const vector& b) const { a -= b; }
};

/** \brief In-place mutate right operand for vector subtraction. */
struct DgCellVecMutSubB
{
    inline void operator()(const vector& a, vector& b) const { b = a - b; }
};

/** \brief In-place mutate left operand for vector cross product. */
struct DgCellVecMutCrossA
{
    inline void operator()(vector& a, const vector& b) const { a = a ^ b; }
};

/** \brief In-place mutate right operand for vector cross product. */
struct DgCellVecMutCrossB
{
    inline void operator()(const vector& a, vector& b) const { b = a ^ b; }
};

/** \brief In-place mutate left operand for `&` with tensor-like types. */
struct DgCellVecMutAndA
{
    template<class T>
    inline void operator()(vector& a, const T& b) const { a = a & b; }
};

/** \brief In-place mutate right operand for `&` with tensor-like types. */
struct DgCellVecMutAndB
{
    template<class T>
    inline void operator()(const T& a, vector& b) const { b = a & b; }
};

/** \brief Binary vector-vector operators: `+`, `-`, `^`, `&`, `*`. */
// Vector - Vector
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalAdd,   DgCellVecMutAddA,   DgCellVecMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalSub,   DgCellVecMutSubA,   DgCellVecMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(^, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalCross, DgCellVecMutCrossA, DgCellVecMutCrossB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<scalar>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalDot,   Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(*, cellGaussField<tensor>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalOuter, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

/** \brief Unary vector operators: `mag`, `magSqr`. */
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(mag,    cellGaussField<scalar>, cellGaussField<vector>, DgCellVecEvalMag,    Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(magSqr, cellGaussField<scalar>, cellGaussField<vector>, DgCellVecEvalMagSqr, Foam::dgCellTmp::NoMutateUnary)

/** \brief Binary operators for `vector & tensor`. */
// Vector - Tensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<tensor>, DgCellVecEvalAndVT, DgCellVecMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<tensor>, cellGaussField<vector>, DgCellVecEvalAndTV, Foam::dgCellTmp::NoMutateLeft, DgCellVecMutAndB)

/** \brief Binary operators for `vector & symmTensor`. */
// Vector - symmTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<symmTensor>, DgCellVecEvalAndVT, DgCellVecMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<symmTensor>, cellGaussField<vector>, DgCellVecEvalAndTV, Foam::dgCellTmp::NoMutateLeft, DgCellVecMutAndB)

/** \brief Binary operators for `vector & sphericalTensor`. */
// Vector - sphericalTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<sphericalTensor>, DgCellVecEvalAndVT, DgCellVecMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<sphericalTensor>, cellGaussField<vector>, DgCellVecEvalAndTV, Foam::dgCellTmp::NoMutateLeft, DgCellVecMutAndB)

} // End namespace Foam

#endif

// ************************************************************************* //
