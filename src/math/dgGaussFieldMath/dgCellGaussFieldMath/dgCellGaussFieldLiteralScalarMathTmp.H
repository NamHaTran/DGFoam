/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Operators between literal scalar and cellGaussField using tmp-return API.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellGaussFieldLiteralScalarMathTmp_H
#define Foam_dgCellGaussFieldLiteralScalarMathTmp_H

#include "cellGaussField.H"
#include "dgCellGaussFieldTmpOpCore.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include <cmath>

namespace Foam
{

template<class Type>
struct DgMulFieldScalarEval
{
    inline Type operator()(const Type& a, const scalar& b) const { return a*b; }
};

template<class Type>
struct DgMulScalarFieldEval
{
    inline Type operator()(const scalar& a, const Type& b) const { return a*b; }
};

template<class Type>
struct DgMulInPlace
{
    inline void operator()(Type& a, const scalar& b) const { a *= b; }
    inline void operator()(const scalar& a, Type& b) const { b *= a; }
};

template<class Type>
struct DgDivFieldScalarEval
{
    inline Type operator()(const Type& a, const scalar& b) const { return a/b; }
};

template<class Type>
struct DgDivFieldScalarInPlace
{
    inline void operator()(Type& a, const scalar& b) const { a /= b; }
};

struct DgDivScalarFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a/b; }
};

struct DgDivScalarFieldScalarInPlace
{
    inline void operator()(const scalar& a, scalar& b) const { b = a/b; }
};

struct DgAddFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a + b; }
};

struct DgAddFieldScalarInPlace
{
    inline void operator()(scalar& a, const scalar& b) const { a += b; }
};

struct DgAddScalarFieldEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a + b; }
};

struct DgAddScalarFieldInPlace
{
    inline void operator()(const scalar& a, scalar& b) const { b = a + b; }
};

struct DgSubFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a - b; }
};

struct DgSubFieldScalarInPlace
{
    inline void operator()(scalar& a, const scalar& b) const { a -= b; }
};

struct DgSubScalarFieldEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a - b; }
};

struct DgSubScalarFieldInPlace
{
    inline void operator()(const scalar& a, scalar& b) const { b = a - b; }
};

struct DgPowEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const
    {
        return std::pow(a, b);
    }
};

struct DgPowMutate
{
    inline void operator()(scalar& a, const scalar& b) const
    {
        a = std::pow(a, b);
    }
};

struct DgSqrtEval
{
    inline scalar operator()(const scalar& a) const
    {
        return std::sqrt(a);
    }
};

struct DgSqrtMutate
{
    inline void operator()(scalar& a) const
    {
        a = std::sqrt(a);
    }
};

struct DgLessPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

struct DgLessEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

struct DgGreaterPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

struct DgGreaterEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

struct DgEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

struct DgNePred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

// scalar * field and field * scalar
DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellGaussField<scalar>,
    scalar,
    cellGaussField<scalar>,
    DgMulScalarFieldEval<scalar>,
    DgMulInPlace<scalar>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    scalar,
    DgMulFieldScalarEval<scalar>,
    DgMulInPlace<scalar>
)

DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellGaussField<vector>,
    scalar,
    cellGaussField<vector>,
    DgMulScalarFieldEval<vector>,
    DgMulInPlace<vector>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellGaussField<vector>,
    cellGaussField<vector>,
    scalar,
    DgMulFieldScalarEval<vector>,
    DgMulInPlace<vector>
)

DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellGaussField<tensor>,
    scalar,
    cellGaussField<tensor>,
    DgMulScalarFieldEval<tensor>,
    DgMulInPlace<tensor>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellGaussField<tensor>,
    cellGaussField<tensor>,
    scalar,
    DgMulFieldScalarEval<tensor>,
    DgMulInPlace<tensor>
)

DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellGaussField<symmTensor>,
    scalar,
    cellGaussField<symmTensor>,
    DgMulScalarFieldEval<symmTensor>,
    DgMulInPlace<symmTensor>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellGaussField<symmTensor>,
    cellGaussField<symmTensor>,
    scalar,
    DgMulFieldScalarEval<symmTensor>,
    DgMulInPlace<symmTensor>
)

DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    cellGaussField<sphericalTensor>,
    scalar,
    cellGaussField<sphericalTensor>,
    DgMulScalarFieldEval<sphericalTensor>,
    DgMulInPlace<sphericalTensor>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    cellGaussField<sphericalTensor>,
    cellGaussField<sphericalTensor>,
    scalar,
    DgMulFieldScalarEval<sphericalTensor>,
    DgMulInPlace<sphericalTensor>
)

// division
DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    /,
    cellGaussField<scalar>,
    scalar,
    cellGaussField<scalar>,
    DgDivScalarFieldScalarEval,
    DgDivScalarFieldScalarInPlace
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    scalar,
    DgDivFieldScalarEval<scalar>,
    DgDivFieldScalarInPlace<scalar>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellGaussField<vector>,
    cellGaussField<vector>,
    scalar,
    DgDivFieldScalarEval<vector>,
    DgDivFieldScalarInPlace<vector>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellGaussField<tensor>,
    cellGaussField<tensor>,
    scalar,
    DgDivFieldScalarEval<tensor>,
    DgDivFieldScalarInPlace<tensor>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellGaussField<symmTensor>,
    cellGaussField<symmTensor>,
    scalar,
    DgDivFieldScalarEval<symmTensor>,
    DgDivFieldScalarInPlace<symmTensor>
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    cellGaussField<sphericalTensor>,
    cellGaussField<sphericalTensor>,
    scalar,
    DgDivFieldScalarEval<sphericalTensor>,
    DgDivFieldScalarInPlace<sphericalTensor>
)

// +/- for scalar fields
DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    +,
    cellGaussField<scalar>,
    scalar,
    cellGaussField<scalar>,
    DgAddScalarFieldEval,
    DgAddScalarFieldInPlace
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    +,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    scalar,
    DgAddFieldScalarEval,
    DgAddFieldScalarInPlace
)

DG_CELL_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    -,
    cellGaussField<scalar>,
    scalar,
    cellGaussField<scalar>,
    DgSubScalarFieldEval,
    DgSubScalarFieldInPlace
)

DG_CELL_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    -,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    scalar,
    DgSubFieldScalarEval,
    DgSubFieldScalarInPlace
)

// pow(field<scalar>, scalar)
DG_CELL_TMP_DEFINE_FIELD_SCALAR_FUNCTION
(
    pow,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    scalar,
    DgPowEval,
    DgPowMutate
)

// sqrt(field<scalar>)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION
(
    sqrt,
    cellGaussField<scalar>,
    cellGaussField<scalar>,
    DgSqrtEval,
    DgSqrtMutate
)

// comparisons: field-scalar and scalar-field
DG_CELL_TMP_DEFINE_COMPARE_FIELD_SCALAR(<,  cellGaussField<scalar>, scalar, DgLessPred)
DG_CELL_TMP_DEFINE_COMPARE_FIELD_SCALAR(<=, cellGaussField<scalar>, scalar, DgLessEqPred)
DG_CELL_TMP_DEFINE_COMPARE_FIELD_SCALAR(>,  cellGaussField<scalar>, scalar, DgGreaterPred)
DG_CELL_TMP_DEFINE_COMPARE_FIELD_SCALAR(>=, cellGaussField<scalar>, scalar, DgGreaterEqPred)
DG_CELL_TMP_DEFINE_COMPARE_FIELD_SCALAR(==, cellGaussField<scalar>, scalar, DgEqPred)
DG_CELL_TMP_DEFINE_COMPARE_FIELD_SCALAR(!=, cellGaussField<scalar>, scalar, DgNePred)

DG_CELL_TMP_DEFINE_COMPARE_SCALAR_FIELD(<,  scalar, cellGaussField<scalar>, DgLessPred)
DG_CELL_TMP_DEFINE_COMPARE_SCALAR_FIELD(<=, scalar, cellGaussField<scalar>, DgLessEqPred)
DG_CELL_TMP_DEFINE_COMPARE_SCALAR_FIELD(>,  scalar, cellGaussField<scalar>, DgGreaterPred)
DG_CELL_TMP_DEFINE_COMPARE_SCALAR_FIELD(>=, scalar, cellGaussField<scalar>, DgGreaterEqPred)
DG_CELL_TMP_DEFINE_COMPARE_SCALAR_FIELD(==, scalar, cellGaussField<scalar>, DgEqPred)
DG_CELL_TMP_DEFINE_COMPARE_SCALAR_FIELD(!=, scalar, cellGaussField<scalar>, DgNePred)

} // End namespace Foam

#endif

// ************************************************************************* //
