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
    Vector math operations between cellDof objects (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofVectorMathTmp_H
#define Foam_dgCellDofVectorMathTmp_H

#include "cellDof.H"
#include "VectorSpace.H"
#include "dgCellDofTmpOpCore.H"

namespace Foam
{

struct DgCellDofVecEvalAdd
{
    inline vector operator()(const vector& a, const vector& b) const { return a + b; }
};

struct DgCellDofVecEvalSub
{
    inline vector operator()(const vector& a, const vector& b) const { return a - b; }
};

struct DgCellDofVecEvalCross
{
    inline vector operator()(const vector& a, const vector& b) const { return a ^ b; }
};

struct DgCellDofVecEvalDot
{
    inline scalar operator()(const vector& a, const vector& b) const { return a & b; }
};

struct DgCellDofVecEvalOuter
{
    inline tensor operator()(const vector& a, const vector& b) const { return a * b; }
};

struct DgCellDofVecEvalAndVT
{
    template<class T>
    inline vector operator()(const vector& a, const T& b) const { return a & b; }
};

struct DgCellDofVecEvalAndTV
{
    template<class T>
    inline vector operator()(const T& a, const vector& b) const { return a & b; }
};

struct DgCellDofVecEvalMag
{
    inline scalar operator()(const vector& a) const { return mag(a); }
};

struct DgCellDofVecEvalMagSqr
{
    inline scalar operator()(const vector& a) const { return magSqr(a); }
};

struct DgCellDofVecMutAddA
{
    inline void operator()(vector& a, const vector& b) const { a += b; }
};

struct DgCellDofVecMutAddB
{
    inline void operator()(const vector& a, vector& b) const { b += a; }
};

struct DgCellDofVecMutSubA
{
    inline void operator()(vector& a, const vector& b) const { a -= b; }
};

struct DgCellDofVecMutSubB
{
    inline void operator()(const vector& a, vector& b) const { b = a - b; }
};

struct DgCellDofVecMutCrossA
{
    inline void operator()(vector& a, const vector& b) const { a = a ^ b; }
};

struct DgCellDofVecMutCrossB
{
    inline void operator()(const vector& a, vector& b) const { b = a ^ b; }
};

struct DgCellDofVecMutAndA
{
    template<class T>
    inline void operator()(vector& a, const T& b) const { a = a & b; }
};

struct DgCellDofVecMutAndB
{
    template<class T>
    inline void operator()(const T& a, vector& b) const { b = a & b; }
};

// Vector - Vector
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(+, cellDof<vector>, cellDof<vector>, cellDof<vector>, DgCellDofVecEvalAdd,   DgCellDofVecMutAddA,   DgCellDofVecMutAddB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(-, cellDof<vector>, cellDof<vector>, cellDof<vector>, DgCellDofVecEvalSub,   DgCellDofVecMutSubA,   DgCellDofVecMutSubB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(^, cellDof<vector>, cellDof<vector>, cellDof<vector>, DgCellDofVecEvalCross, DgCellDofVecMutCrossA, DgCellDofVecMutCrossB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&, cellDof<scalar>, cellDof<vector>, cellDof<vector>, DgCellDofVecEvalDot,   Foam::dgCellDofTmp::NoMutateLeft, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(*, cellDof<tensor>, cellDof<vector>, cellDof<vector>, DgCellDofVecEvalOuter, Foam::dgCellDofTmp::NoMutateLeft, Foam::dgCellDofTmp::NoMutateRight)

// Unary vector ops
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(mag,    cellDof<scalar>, cellDof<vector>, DgCellDofVecEvalMag,    Foam::dgCellDofTmp::NoMutateUnary)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(magSqr, cellDof<scalar>, cellDof<vector>, DgCellDofVecEvalMagSqr, Foam::dgCellDofTmp::NoMutateUnary)

// Vector - Tensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&, cellDof<vector>, cellDof<vector>, cellDof<tensor>, DgCellDofVecEvalAndVT, DgCellDofVecMutAndA, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&, cellDof<vector>, cellDof<tensor>, cellDof<vector>, DgCellDofVecEvalAndTV, Foam::dgCellDofTmp::NoMutateLeft, DgCellDofVecMutAndB)

// Vector - symmTensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&, cellDof<vector>, cellDof<vector>, cellDof<symmTensor>, DgCellDofVecEvalAndVT, DgCellDofVecMutAndA, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&, cellDof<vector>, cellDof<symmTensor>, cellDof<vector>, DgCellDofVecEvalAndTV, Foam::dgCellDofTmp::NoMutateLeft, DgCellDofVecMutAndB)

// Vector - sphericalTensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&, cellDof<vector>, cellDof<vector>, cellDof<sphericalTensor>, DgCellDofVecEvalAndVT, DgCellDofVecMutAndA, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&, cellDof<vector>, cellDof<sphericalTensor>, cellDof<vector>, DgCellDofVecEvalAndTV, Foam::dgCellDofTmp::NoMutateLeft, DgCellDofVecMutAndB)

} // End namespace Foam

#endif

// ************************************************************************* //

