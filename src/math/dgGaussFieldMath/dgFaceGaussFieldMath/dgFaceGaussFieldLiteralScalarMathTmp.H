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
    Operators between literal scalar and faceGaussField using tmp-return API.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFaceGaussFieldLiteralScalarMathTmp_H
#define Foam_dgFaceGaussFieldLiteralScalarMathTmp_H

#include "faceGaussField.H"
#include "dgFaceGaussFieldTmpOpCore.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include <cmath>

namespace Foam
{

/** \brief Evaluate `field * scalar` at one face Gauss point. */
template<class Type>
struct DgFaceMulFieldScalarEval
{
    inline Type operator()(const Type& a, const scalar& b) const { return a*b; }
};

/** \brief Evaluate `scalar * field` at one face Gauss point. */
template<class Type>
struct DgFaceMulScalarFieldEval
{
    inline Type operator()(const scalar& a, const Type& b) const { return a*b; }
};

/** \brief In-place mutator for scalar multiplication on face fields. */
template<class Type>
struct DgFaceMulInPlace
{
    inline void operator()(Type& ap, Type& am, const scalar& b) const
    {
        ap *= b;
        am *= b;
    }

    inline void operator()(const scalar& a, Type& bp, Type& bm) const
    {
        bp *= a;
        bm *= a;
    }
};

/** \brief Evaluate `field / scalar` at one face Gauss point. */
template<class Type>
struct DgFaceDivFieldScalarEval
{
    inline Type operator()(const Type& a, const scalar& b) const { return a/b; }
};

/** \brief In-place mutator for `field / scalar` on face fields. */
template<class Type>
struct DgFaceDivFieldScalarInPlace
{
    inline void operator()(Type& ap, Type& am, const scalar& b) const
    {
        ap /= b;
        am /= b;
    }
};

/** \brief Evaluate `scalar / field` for scalar face fields. */
struct DgFaceDivScalarFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a/b; }
};

/** \brief In-place mutator for `scalar / scalarFaceField`. */
struct DgFaceDivScalarFieldScalarInPlace
{
    inline void operator()(const scalar& a, scalar& bp, scalar& bm) const
    {
        bp = a/bp;
        bm = a/bm;
    }
};

/** \brief Evaluate `field + scalar` for scalar face fields. */
struct DgFaceAddFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a + b; }
};

/** \brief In-place mutator for `field + scalar`. */
struct DgFaceAddFieldScalarInPlace
{
    inline void operator()(scalar& ap, scalar& am, const scalar& b) const
    {
        ap += b;
        am += b;
    }
};

/** \brief Evaluate `scalar + field` for scalar face fields. */
struct DgFaceAddScalarFieldEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a + b; }
};

/** \brief In-place mutator for `scalar + field`. */
struct DgFaceAddScalarFieldInPlace
{
    inline void operator()(const scalar& a, scalar& bp, scalar& bm) const
    {
        bp = a + bp;
        bm = a + bm;
    }
};

/** \brief Evaluate `field - scalar` for scalar face fields. */
struct DgFaceSubFieldScalarEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a - b; }
};

/** \brief In-place mutator for `field - scalar`. */
struct DgFaceSubFieldScalarInPlace
{
    inline void operator()(scalar& ap, scalar& am, const scalar& b) const
    {
        ap -= b;
        am -= b;
    }
};

/** \brief Evaluate `scalar - field` for scalar face fields. */
struct DgFaceSubScalarFieldEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const { return a - b; }
};

/** \brief In-place mutator for `scalar - field`. */
struct DgFaceSubScalarFieldInPlace
{
    inline void operator()(const scalar& a, scalar& bp, scalar& bm) const
    {
        bp = a - bp;
        bm = a - bm;
    }
};

/** \brief Evaluate `pow(fieldScalar, literalScalar)`. */
struct DgFacePowEval
{
    inline scalar operator()(const scalar& a, const scalar& b) const
    {
        return std::pow(a, b);
    }
};

/** \brief In-place mutator for `pow(fieldScalar, literalScalar)`. */
struct DgFacePowMutate
{
    inline void operator()(scalar& ap, scalar& am, const scalar& b) const
    {
        ap = std::pow(ap, b);
        am = std::pow(am, b);
    }
};

/** \brief Evaluate `sqrt(fieldScalar)` at one face Gauss point. */
struct DgFaceSqrtEval
{
    inline scalar operator()(const scalar& a) const
    {
        return std::sqrt(a);
    }
};

/** \brief In-place mutator for `sqrt(fieldScalar)`. */
struct DgFaceSqrtMutate
{
    inline void operator()(scalar& ap, scalar& am) const
    {
        ap = std::sqrt(ap);
        am = std::sqrt(am);
    }
};

/** \brief Predicate functor for `<` against scalar literal. */
struct DgFaceLessPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a < b; }
};

/** \brief Predicate functor for `<=` against scalar literal. */
struct DgFaceLessEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a <= b; }
};

/** \brief Predicate functor for `>` against scalar literal. */
struct DgFaceGreaterPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a > b; }
};

/** \brief Predicate functor for `>=` against scalar literal. */
struct DgFaceGreaterEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a >= b; }
};

/** \brief Predicate functor for `==` against scalar literal. */
struct DgFaceEqPred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a == b; }
};

/** \brief Predicate functor for `!=` against scalar literal. */
struct DgFaceNePred
{
    inline bool operator()(const scalar& a, const scalar& b) const { return a != b; }
};

/** \brief Multiplication operators between scalar literal and face field types. */
// scalar * field and field * scalar
DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    faceGaussField<scalar>,
    scalar,
    faceGaussField<scalar>,
    DgFaceMulScalarFieldEval<scalar>,
    DgFaceMulInPlace<scalar>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    scalar,
    DgFaceMulFieldScalarEval<scalar>,
    DgFaceMulInPlace<scalar>
)

DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    faceGaussField<vector>,
    scalar,
    faceGaussField<vector>,
    DgFaceMulScalarFieldEval<vector>,
    DgFaceMulInPlace<vector>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    faceGaussField<vector>,
    faceGaussField<vector>,
    scalar,
    DgFaceMulFieldScalarEval<vector>,
    DgFaceMulInPlace<vector>
)

DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    faceGaussField<tensor>,
    scalar,
    faceGaussField<tensor>,
    DgFaceMulScalarFieldEval<tensor>,
    DgFaceMulInPlace<tensor>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    faceGaussField<tensor>,
    faceGaussField<tensor>,
    scalar,
    DgFaceMulFieldScalarEval<tensor>,
    DgFaceMulInPlace<tensor>
)

DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    faceGaussField<symmTensor>,
    scalar,
    faceGaussField<symmTensor>,
    DgFaceMulScalarFieldEval<symmTensor>,
    DgFaceMulInPlace<symmTensor>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    faceGaussField<symmTensor>,
    faceGaussField<symmTensor>,
    scalar,
    DgFaceMulFieldScalarEval<symmTensor>,
    DgFaceMulInPlace<symmTensor>
)

DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    *,
    faceGaussField<sphericalTensor>,
    scalar,
    faceGaussField<sphericalTensor>,
    DgFaceMulScalarFieldEval<sphericalTensor>,
    DgFaceMulInPlace<sphericalTensor>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    *,
    faceGaussField<sphericalTensor>,
    faceGaussField<sphericalTensor>,
    scalar,
    DgFaceMulFieldScalarEval<sphericalTensor>,
    DgFaceMulInPlace<sphericalTensor>
)

// division
DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    /,
    faceGaussField<scalar>,
    scalar,
    faceGaussField<scalar>,
    DgFaceDivScalarFieldScalarEval,
    DgFaceDivScalarFieldScalarInPlace
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    scalar,
    DgFaceDivFieldScalarEval<scalar>,
    DgFaceDivFieldScalarInPlace<scalar>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    faceGaussField<vector>,
    faceGaussField<vector>,
    scalar,
    DgFaceDivFieldScalarEval<vector>,
    DgFaceDivFieldScalarInPlace<vector>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    faceGaussField<tensor>,
    faceGaussField<tensor>,
    scalar,
    DgFaceDivFieldScalarEval<tensor>,
    DgFaceDivFieldScalarInPlace<tensor>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    faceGaussField<symmTensor>,
    faceGaussField<symmTensor>,
    scalar,
    DgFaceDivFieldScalarEval<symmTensor>,
    DgFaceDivFieldScalarInPlace<symmTensor>
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    /,
    faceGaussField<sphericalTensor>,
    faceGaussField<sphericalTensor>,
    scalar,
    DgFaceDivFieldScalarEval<sphericalTensor>,
    DgFaceDivFieldScalarInPlace<sphericalTensor>
)

// +/- for scalar fields
DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    +,
    faceGaussField<scalar>,
    scalar,
    faceGaussField<scalar>,
    DgFaceAddScalarFieldEval,
    DgFaceAddScalarFieldInPlace
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    +,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    scalar,
    DgFaceAddFieldScalarEval,
    DgFaceAddFieldScalarInPlace
)

DG_FACE_TMP_DEFINE_SCALAR_FIELD_OPERATOR
(
    -,
    faceGaussField<scalar>,
    scalar,
    faceGaussField<scalar>,
    DgFaceSubScalarFieldEval,
    DgFaceSubScalarFieldInPlace
)

DG_FACE_TMP_DEFINE_FIELD_SCALAR_OPERATOR
(
    -,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    scalar,
    DgFaceSubFieldScalarEval,
    DgFaceSubFieldScalarInPlace
)

// pow(field<scalar>, scalar)
DG_FACE_TMP_DEFINE_FIELD_SCALAR_FUNCTION
(
    pow,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    scalar,
    DgFacePowEval,
    DgFacePowMutate
)

// sqrt(field<scalar>)
DG_FACE_TMP_DEFINE_UNARY_FUNCTION
(
    sqrt,
    faceGaussField<scalar>,
    faceGaussField<scalar>,
    DgFaceSqrtEval,
    DgFaceSqrtMutate
)

// comparisons: field-scalar and scalar-field
DG_FACE_TMP_DEFINE_COMPARE_FIELD_SCALAR(<,  faceGaussField<scalar>, scalar, DgFaceLessPred)
DG_FACE_TMP_DEFINE_COMPARE_FIELD_SCALAR(<=, faceGaussField<scalar>, scalar, DgFaceLessEqPred)
DG_FACE_TMP_DEFINE_COMPARE_FIELD_SCALAR(>,  faceGaussField<scalar>, scalar, DgFaceGreaterPred)
DG_FACE_TMP_DEFINE_COMPARE_FIELD_SCALAR(>=, faceGaussField<scalar>, scalar, DgFaceGreaterEqPred)
DG_FACE_TMP_DEFINE_COMPARE_FIELD_SCALAR(==, faceGaussField<scalar>, scalar, DgFaceEqPred)
DG_FACE_TMP_DEFINE_COMPARE_FIELD_SCALAR(!=, faceGaussField<scalar>, scalar, DgFaceNePred)

DG_FACE_TMP_DEFINE_COMPARE_SCALAR_FIELD(<,  scalar, faceGaussField<scalar>, DgFaceLessPred)
DG_FACE_TMP_DEFINE_COMPARE_SCALAR_FIELD(<=, scalar, faceGaussField<scalar>, DgFaceLessEqPred)
DG_FACE_TMP_DEFINE_COMPARE_SCALAR_FIELD(>,  scalar, faceGaussField<scalar>, DgFaceGreaterPred)
DG_FACE_TMP_DEFINE_COMPARE_SCALAR_FIELD(>=, scalar, faceGaussField<scalar>, DgFaceGreaterEqPred)
DG_FACE_TMP_DEFINE_COMPARE_SCALAR_FIELD(==, scalar, faceGaussField<scalar>, DgFaceEqPred)
DG_FACE_TMP_DEFINE_COMPARE_SCALAR_FIELD(!=, scalar, faceGaussField<scalar>, DgFaceNePred)

} // End namespace Foam

#endif

// ************************************************************************* //
