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
    Operators between literal scalar and cellDof using tmp-return API.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofLiteralScalarMathTmp_H
#define Foam_dgCellDofLiteralScalarMathTmp_H

#include "cellDof.H"
#include "dgCellDofTmpOpCore.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include <cmath>

namespace Foam
{

template<class Type>
struct DgCellDofMulFieldScalarEval
{
    inline Type operator()(const Type& a, const scalar& b) const { return a*b; }
};

template<class Type>
struct DgCellDofMulScalarFieldEval
{
    inline Type operator()(const scalar& a, const Type& b) const { return a*b; }
};

template<class Type>
struct DgCellDofMulInPlace
{
    inline void operator()(Type& a, const scalar& b) const { a *= b; }
    inline void operator()(const scalar& a, Type& b) const { b *= a; }
};

template<class Type>
struct DgCellDofDivFieldScalarEval
{
    inline Type operator()(const Type& a, const scalar& b) const { return a/b; }
};

template<class Type>
struct DgCellDofDivFieldScalarInPlace
{
    inline void operator()(Type& a, const scalar& b) const { a /= b; }
};

struct DgCellDofDivScalarFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a/b; }
};

struct DgCellDofDivScalarFieldScalarInPlace
{
    inline void operator()(const scalar& a, scalar& b) const { b = a/b; }
};

struct DgCellDofAddFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a + b; }
};

struct DgCellDofAddFieldScalarInPlace
{
    inline void operator()(scalar& a, const scalar& b) const { a += b; }
};

struct DgCellDofAddScalarFieldEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a + b; }
};

struct DgCellDofAddScalarFieldInPlace
{
    inline void operator()(const scalar& a, scalar& b) const { b = a + b; }
};

struct DgCellDofSubFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a - b; }
};

struct DgCellDofSubFieldScalarInPlace
{
    inline void operator()(scalar& a, const scalar& b) const { a -= b; }
};

struct DgCellDofSubScalarFieldEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a - b; }
};

struct DgCellDofSubScalarFieldInPlace
{
    inline void operator()(const scalar& a, scalar& b) const { b = a - b; }
};

struct DgCellDofPowEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const
    {
        return std::pow(a, b);
    }
};

struct DgCellDofPowMutate
{
    inline void operator()(scalar& a, const scalar& b) const
    {
        a = std::pow(a, b);
    }
};

struct DgCellDofSqrtEval
{
    inline scalar operator()(const scalar& a) const
    {
        return std::sqrt(a);
    }
};

struct DgCellDofSqrtMutate
{
    inline void operator()(scalar& a) const
    {
        a = std::sqrt(a);
    }
};

struct DgCellDofLessPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

struct DgCellDofLessEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

struct DgCellDofGreaterPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

struct DgCellDofGreaterEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

struct DgCellDofEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

struct DgCellDofNePred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

// scalar * field and field * scalar
DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellDof<scalar>,
    scalar,
    cellDof<scalar>,
    DgCellDofMulScalarFieldEval<scalar>,
    DgCellDofMulInPlace<scalar>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellDof<scalar>,
    cellDof<scalar>,
    scalar,
    DgCellDofMulFieldScalarEval<scalar>,
    DgCellDofMulInPlace<scalar>
)

DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellDof<vector>,
    scalar,
    cellDof<vector>,
    DgCellDofMulScalarFieldEval<vector>,
    DgCellDofMulInPlace<vector>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellDof<vector>,
    cellDof<vector>,
    scalar,
    DgCellDofMulFieldScalarEval<vector>,
    DgCellDofMulInPlace<vector>
)

DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellDof<tensor>,
    scalar,
    cellDof<tensor>,
    DgCellDofMulScalarFieldEval<tensor>,
    DgCellDofMulInPlace<tensor>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellDof<tensor>,
    cellDof<tensor>,
    scalar,
    DgCellDofMulFieldScalarEval<tensor>,
    DgCellDofMulInPlace<tensor>
)

DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellDof<symmTensor>,
    scalar,
    cellDof<symmTensor>,
    DgCellDofMulScalarFieldEval<symmTensor>,
    DgCellDofMulInPlace<symmTensor>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellDof<symmTensor>,
    cellDof<symmTensor>,
    scalar,
    DgCellDofMulFieldScalarEval<symmTensor>,
    DgCellDofMulInPlace<symmTensor>
)

DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellDof<sphericalTensor>,
    scalar,
    cellDof<sphericalTensor>,
    DgCellDofMulScalarFieldEval<sphericalTensor>,
    DgCellDofMulInPlace<sphericalTensor>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellDof<sphericalTensor>,
    cellDof<sphericalTensor>,
    scalar,
    DgCellDofMulFieldScalarEval<sphericalTensor>,
    DgCellDofMulInPlace<sphericalTensor>
)

// division
DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    /,
    cellDof<scalar>,
    scalar,
    cellDof<scalar>,
    DgCellDofDivScalarFieldScalarEval,
    DgCellDofDivScalarFieldScalarInPlace
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellDof<scalar>,
    cellDof<scalar>,
    scalar,
    DgCellDofDivFieldScalarEval<scalar>,
    DgCellDofDivFieldScalarInPlace<scalar>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellDof<vector>,
    cellDof<vector>,
    scalar,
    DgCellDofDivFieldScalarEval<vector>,
    DgCellDofDivFieldScalarInPlace<vector>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellDof<tensor>,
    cellDof<tensor>,
    scalar,
    DgCellDofDivFieldScalarEval<tensor>,
    DgCellDofDivFieldScalarInPlace<tensor>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellDof<symmTensor>,
    cellDof<symmTensor>,
    scalar,
    DgCellDofDivFieldScalarEval<symmTensor>,
    DgCellDofDivFieldScalarInPlace<symmTensor>
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellDof<sphericalTensor>,
    cellDof<sphericalTensor>,
    scalar,
    DgCellDofDivFieldScalarEval<sphericalTensor>,
    DgCellDofDivFieldScalarInPlace<sphericalTensor>
)

// +/- for scalar fields
DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    +,
    cellDof<scalar>,
    scalar,
    cellDof<scalar>,
    DgCellDofAddScalarFieldEval,
    DgCellDofAddScalarFieldInPlace
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    +,
    cellDof<scalar>,
    cellDof<scalar>,
    scalar,
    DgCellDofAddFieldScalarEval,
    DgCellDofAddFieldScalarInPlace
)

DG_CELL_DOF_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    -,
    cellDof<scalar>,
    scalar,
    cellDof<scalar>,
    DgCellDofSubScalarFieldEval,
    DgCellDofSubScalarFieldInPlace
)

DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    -,
    cellDof<scalar>,
    cellDof<scalar>,
    scalar,
    DgCellDofSubFieldScalarEval,
    DgCellDofSubFieldScalarInPlace
)

// pow(field<scalar>, scalar)
DG_CELL_DOF_TMP_DEFINE_FIELD_SCALAR_FUNCTION
(
    pow,
    cellDof<scalar>,
    cellDof<scalar>,
    scalar,
    DgCellDofPowEval,
    DgCellDofPowMutate
)

// sqrt(field<scalar>)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION
(
    sqrt,
    cellDof<scalar>,
    cellDof<scalar>,
    DgCellDofSqrtEval,
    DgCellDofSqrtMutate
)

// comparisons: field-scalar
DG_CELL_DOF_TMP_DEFINE_COMPARE_FIELD_SCALAR(<,  cellDof<scalar>, scalar, DgCellDofLessPred)
DG_CELL_DOF_TMP_DEFINE_COMPARE_FIELD_SCALAR(<=, cellDof<scalar>, scalar, DgCellDofLessEqPred)
DG_CELL_DOF_TMP_DEFINE_COMPARE_FIELD_SCALAR(>,  cellDof<scalar>, scalar, DgCellDofGreaterPred)
DG_CELL_DOF_TMP_DEFINE_COMPARE_FIELD_SCALAR(>=, cellDof<scalar>, scalar, DgCellDofGreaterEqPred)
DG_CELL_DOF_TMP_DEFINE_COMPARE_FIELD_SCALAR(==, cellDof<scalar>, scalar, DgCellDofEqPred)
DG_CELL_DOF_TMP_DEFINE_COMPARE_FIELD_SCALAR(!=, cellDof<scalar>, scalar, DgCellDofNePred)

} // End namespace Foam

#endif

// ************************************************************************* //

