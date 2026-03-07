/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

struct DgCellEvalDev
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(dev(a)) { return dev(a); }
};

struct DgCellEvalSymm
{
    inline symmTensor operator()(const tensor& a) const { return symm(a); }
};

struct DgCellEvalSkew
{
    inline tensor operator()(const tensor& a) const { return skew(a); }
};

struct DgCellEvalTr
{
    template<class A>
    inline scalar operator()(const A& a) const { return tr(a); }
};

struct DgCellEvalDet
{
    inline scalar operator()(const tensor& a) const { return det(a); }
};

struct DgCellEvalT
{
    inline tensor operator()(const tensor& a) const { return T(a); }
};

struct DgCellEvalIdentitySymm
{
    inline symmTensor operator()(const symmTensor& a) const { return a; }
};

struct DgCellEvalInv
{
    template<class A>
    inline tensor operator()(const A& a) const { return inv(a); }
};

struct DgCellMutDevTensor
{
    inline void operator()(tensor& a) const { a = dev(a); }
};

struct DgCellMutSkewTensor
{
    inline void operator()(tensor& a) const { a = skew(a); }
};

struct DgCellMutTTensor
{
    inline void operator()(tensor& a) const { a = T(a); }
};

struct DgCellMutDevSymm
{
    inline void operator()(symmTensor& a) const { a = dev(a); }
};

struct DgCellMutNoopSymm
{
    inline void operator()(symmTensor&) const {}
};

struct DgCellEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

struct DgCellEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

struct DgCellEvalAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b)
    {
        return a & b;
    }
};

struct DgCellEvalAndAnd
{
    template<class A, class B>
    inline scalar operator()(const A& a, const B& b) const
    {
        return a && b;
    }
};

struct DgCellMutAddA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a + b; }
};

struct DgCellMutAddB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a + b; }
};

struct DgCellMutSubA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a - b; }
};

struct DgCellMutSubB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

struct DgCellMutAndA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a & b; }
};

struct DgCellMutAndB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a & b; }
};

// Tensor utilities
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(dev,  cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalDev,       DgCellMutDevTensor)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(symm, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalSymm,      Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(skew, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalSkew,      DgCellMutSkewTensor)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(tr,   cellGaussField<scalar>,     cellGaussField<tensor>,     DgCellEvalTr,        Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(det,  cellGaussField<scalar>,     cellGaussField<tensor>,     DgCellEvalDet,       Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(T,    cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalT,         DgCellMutTTensor)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(inv,  cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalInv,       Foam::dgCellTmp::NoMutateUnary)

// symmTensor utilities
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(tr,   cellGaussField<scalar>,     cellGaussField<symmTensor>, DgCellEvalTr,        Foam::dgCellTmp::NoMutateUnary)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(dev,  cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalDev,       DgCellMutDevSymm)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(T,    cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalIdentitySymm, DgCellMutNoopSymm)
DG_CELL_TMP_DEFINE_UNARY_FUNCTION(inv,  cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalInv,       Foam::dgCellTmp::NoMutateUnary)

// tensor - tensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalAdd,    DgCellMutAddA, DgCellMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalSub,    DgCellMutSubA, DgCellMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalAnd,    DgCellMutAndA, DgCellMutAndB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>, cellGaussField<tensor>,     cellGaussField<tensor>,     DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

// tensor - symmTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalAdd,    DgCellMutAddA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalSub,    DgCellMutSubA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalAnd,    DgCellMutAndA, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>, cellGaussField<tensor>,     cellGaussField<symmTensor>, DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

// symmTensor - tensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<tensor>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalAdd,    Foam::dgCellTmp::NoMutateLeft, DgCellMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<tensor>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalSub,    Foam::dgCellTmp::NoMutateLeft, DgCellMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalAnd,    Foam::dgCellTmp::NoMutateLeft, DgCellMutAndB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>, cellGaussField<symmTensor>, cellGaussField<tensor>,     DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

// symmTensor - symmTensor
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(+,  cellGaussField<symmTensor>, cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalAdd,    DgCellMutAddA, DgCellMutAddB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(-,  cellGaussField<symmTensor>, cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalSub,    DgCellMutSubA, DgCellMutSubB)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&,  cellGaussField<tensor>,     cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalAnd,    Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)
DG_CELL_TMP_DEFINE_BINARY_OPERATOR(&&, cellGaussField<scalar>,     cellGaussField<symmTensor>, cellGaussField<symmTensor>, DgCellEvalAndAnd, Foam::dgCellTmp::NoMutateLeft, Foam::dgCellTmp::NoMutateRight)

} // End namespace Foam

#endif

// ************************************************************************* //

