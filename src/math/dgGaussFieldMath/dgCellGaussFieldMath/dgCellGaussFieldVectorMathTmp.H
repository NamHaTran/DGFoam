/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

struct DgCellVecEvalAdd
{
    inline vector operator()(const vector& a, const vector& b) const { return a + b; }
};

struct DgCellVecEvalSub
{
    inline vector operator()(const vector& a, const vector& b) const { return a - b; }
};

struct DgCellVecEvalCross
{
    inline vector operator()(const vector& a, const vector& b) const { return a ^ b; }
};

struct DgCellVecEvalDot
{
    inline scalar operator()(const vector& a, const vector& b) const { return a & b; }
};

struct DgCellVecEvalOuter
{
    inline tensor operator()(const vector& a, const vector& b) const { return a * b; }
};

struct DgCellVecEvalAndVT
{
    template<class T>
    inline vector operator()(const vector& a, const T& b) const { return a & b; }
};

struct DgCellVecEvalAndTV
{
    template<class T>
    inline vector operator()(const T& a, const vector& b) const { return a & b; }
};

struct DgCellVecEvalMag
{
    inline scalar operator()(const vector& a) const { return mag(a); }
};

struct DgCellVecEvalMagSqr
{
    inline scalar operator()(const vector& a) const { return magSqr(a); }
};

struct DgCellVecMutAddA
{
    inline void operator()(vector& a, const vector& b) const { a += b; }
};

struct DgCellVecMutAddB
{
    inline void operator()(const vector& a, vector& b) const { b += a; }
};

struct DgCellVecMutSubA
{
    inline void operator()(vector& a, const vector& b) const { a -= b; }
};

struct DgCellVecMutSubB
{
    inline void operator()(const vector& a, vector& b) const { b = a - b; }
};

struct DgCellVecMutCrossA
{
    inline void operator()(vector& a, const vector& b) const { a = a ^ b; }
};

struct DgCellVecMutCrossB
{
    inline void operator()(const vector& a, vector& b) const { b = a ^ b; }
};

struct DgCellVecMutAndA
{
    template<class T>
    inline void operator()(vector& a, const T& b) const { a = a & b; }
};

struct DgCellVecMutAndB
{
    template<class T>
    inline void operator()(const T& a, vector& b) const { b = a & b; }
};

// Vector - Vector
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalAdd,   DgCellVecMutAddA,   DgCellVecMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalSub,   DgCellVecMutSubA,   DgCellVecMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(^, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalCross, DgCellVecMutCrossA, DgCellVecMutCrossB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<scalar>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalDot,   Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(*, cellGaussField<tensor>, cellGaussField<vector>, cellGaussField<vector>, DgCellVecEvalOuter, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

DG_CELL_TMP_DEFINE_UNARY_FUNCTION(mag,    cellGaussField<scalar>, cellGaussField<vector>, DgCellVecEvalMag,    Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(magSqr, cellGaussField<scalar>, cellGaussField<vector>, DgCellVecEvalMagSqr, Foam::dgCellTmp::NoMutateUnary)

// Vector - Tensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<tensor>, DgCellVecEvalAndVT, DgCellVecMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<tensor>, cellGaussField<vector>, DgCellVecEvalAndTV, Foam::dgCellTmp::NoMutateLeft, DgCellVecMutAndB)

// Vector - symmTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<symmTensor>, DgCellVecEvalAndVT, DgCellVecMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<symmTensor>, cellGaussField<vector>, DgCellVecEvalAndTV, Foam::dgCellTmp::NoMutateLeft, DgCellVecMutAndB)

// Vector - sphericalTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<vector>, cellGaussField<sphericalTensor>, DgCellVecEvalAndVT, DgCellVecMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&, cellGaussField<vector>, cellGaussField<sphericalTensor>, cellGaussField<vector>, DgCellVecEvalAndTV, Foam::dgCellTmp::NoMutateLeft, DgCellVecMutAndB)

} // End namespace Foam

#endif

// ************************************************************************* //

