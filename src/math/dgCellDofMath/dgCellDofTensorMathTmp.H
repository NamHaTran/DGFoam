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
    Tensor operations between cellDof objects (tmp-aware).

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCellDofTensorMathTmp_H
#define Foam_dgCellDofTensorMathTmp_H

#include "cellDof.H"
#include "Tensor.H"
#include "dgCellDofTmpOpCore.H"

namespace Foam
{

struct DgCellDofTenEvalDev
{
    template<class A>
    inline auto operator()(const A& a) const -> decltype(dev(a)) { return dev(a); }
};

struct DgCellDofTenEvalSymm
{
    inline symmTensor operator()(const tensor& a) const { return symm(a); }
};

struct DgCellDofTenEvalSkew
{
    inline tensor operator()(const tensor& a) const { return skew(a); }
};

struct DgCellDofTenEvalTr
{
    template<class A>
    inline scalar operator()(const A& a) const { return tr(a); }
};

struct DgCellDofTenEvalDet
{
    inline scalar operator()(const tensor& a) const { return det(a); }
};

struct DgCellDofTenEvalT
{
    inline tensor operator()(const tensor& a) const { return T(a); }
};

struct DgCellDofTenEvalIdentitySymm
{
    inline symmTensor operator()(const symmTensor& a) const { return a; }
};

struct DgCellDofTenEvalInv
{
    template<class A>
    inline tensor operator()(const A& a) const { return inv(a); }
};

struct DgCellDofTenMutDevTensor
{
    inline void operator()(tensor& a) const { a = dev(a); }
};

struct DgCellDofTenMutSkewTensor
{
    inline void operator()(tensor& a) const { a = skew(a); }
};

struct DgCellDofTenMutTTensor
{
    inline void operator()(tensor& a) const { a = T(a); }
};

struct DgCellDofTenMutDevSymm
{
    inline void operator()(symmTensor& a) const { a = dev(a); }
};

struct DgCellDofTenMutNoopSymm
{
    inline void operator()(symmTensor&) const {}
};

struct DgCellDofTenEvalAdd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a + b)
    {
        return a + b;
    }
};

struct DgCellDofTenEvalSub
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a - b)
    {
        return a - b;
    }
};

struct DgCellDofTenEvalAnd
{
    template<class A, class B>
    inline auto operator()(const A& a, const B& b) const -> decltype(a & b)
    {
        return a & b;
    }
};

struct DgCellDofTenEvalAndAnd
{
    template<class A, class B>
    inline scalar operator()(const A& a, const B& b) const
    {
        return a && b;
    }
};

struct DgCellDofTenMutAddA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a + b; }
};

struct DgCellDofTenMutAddB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a + b; }
};

struct DgCellDofTenMutSubA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a - b; }
};

struct DgCellDofTenMutSubB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a - b; }
};

struct DgCellDofTenMutAndA
{
    template<class A, class B>
    inline void operator()(A& a, const B& b) const { a = a & b; }
};

struct DgCellDofTenMutAndB
{
    template<class A, class B>
    inline void operator()(const A& a, B& b) const { b = a & b; }
};

// Tensor utilities
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(dev,  cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalDev,       DgCellDofTenMutDevTensor)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(symm, cellDof<symmTensor>, cellDof<tensor>,     DgCellDofTenEvalSymm,      Foam::dgCellDofTmp::NoMutateUnary)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(skew, cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalSkew,      DgCellDofTenMutSkewTensor)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(tr,   cellDof<scalar>,     cellDof<tensor>,     DgCellDofTenEvalTr,        Foam::dgCellDofTmp::NoMutateUnary)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(det,  cellDof<scalar>,     cellDof<tensor>,     DgCellDofTenEvalDet,       Foam::dgCellDofTmp::NoMutateUnary)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(T,    cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalT,         DgCellDofTenMutTTensor)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(inv,  cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalInv,       Foam::dgCellDofTmp::NoMutateUnary)

// symmTensor utilities
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(tr,   cellDof<scalar>,     cellDof<symmTensor>, DgCellDofTenEvalTr,        Foam::dgCellDofTmp::NoMutateUnary)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(dev,  cellDof<symmTensor>, cellDof<symmTensor>, DgCellDofTenEvalDev,       DgCellDofTenMutDevSymm)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(T,    cellDof<symmTensor>, cellDof<symmTensor>, DgCellDofTenEvalIdentitySymm, DgCellDofTenMutNoopSymm)
DG_CELL_DOF_TMP_DEFINE_UNARY_FUNCTION(inv,  cellDof<tensor>,     cellDof<symmTensor>, DgCellDofTenEvalInv,       Foam::dgCellDofTmp::NoMutateUnary)

// tensor - tensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(+,  cellDof<tensor>, cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalAdd,    DgCellDofTenMutAddA, DgCellDofTenMutAddB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(-,  cellDof<tensor>, cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalSub,    DgCellDofTenMutSubA, DgCellDofTenMutSubB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&,  cellDof<tensor>, cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalAnd,    DgCellDofTenMutAndA, DgCellDofTenMutAndB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&&, cellDof<scalar>, cellDof<tensor>,     cellDof<tensor>,     DgCellDofTenEvalAndAnd, Foam::dgCellDofTmp::NoMutateLeft, Foam::dgCellDofTmp::NoMutateRight)

// tensor - symmTensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(+,  cellDof<tensor>, cellDof<tensor>,     cellDof<symmTensor>, DgCellDofTenEvalAdd,    DgCellDofTenMutAddA, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(-,  cellDof<tensor>, cellDof<tensor>,     cellDof<symmTensor>, DgCellDofTenEvalSub,    DgCellDofTenMutSubA, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&,  cellDof<tensor>, cellDof<tensor>,     cellDof<symmTensor>, DgCellDofTenEvalAnd,    DgCellDofTenMutAndA, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&&, cellDof<scalar>, cellDof<tensor>,     cellDof<symmTensor>, DgCellDofTenEvalAndAnd, Foam::dgCellDofTmp::NoMutateLeft, Foam::dgCellDofTmp::NoMutateRight)

// symmTensor - tensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(+,  cellDof<tensor>, cellDof<symmTensor>, cellDof<tensor>,     DgCellDofTenEvalAdd,    Foam::dgCellDofTmp::NoMutateLeft, DgCellDofTenMutAddB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(-,  cellDof<tensor>, cellDof<symmTensor>, cellDof<tensor>,     DgCellDofTenEvalSub,    Foam::dgCellDofTmp::NoMutateLeft, DgCellDofTenMutSubB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&,  cellDof<tensor>, cellDof<symmTensor>, cellDof<tensor>,     DgCellDofTenEvalAnd,    Foam::dgCellDofTmp::NoMutateLeft, DgCellDofTenMutAndB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&&, cellDof<scalar>, cellDof<symmTensor>, cellDof<tensor>,     DgCellDofTenEvalAndAnd, Foam::dgCellDofTmp::NoMutateLeft, Foam::dgCellDofTmp::NoMutateRight)

// symmTensor - symmTensor
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(+,  cellDof<symmTensor>, cellDof<symmTensor>, cellDof<symmTensor>, DgCellDofTenEvalAdd,    DgCellDofTenMutAddA, DgCellDofTenMutAddB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(-,  cellDof<symmTensor>, cellDof<symmTensor>, cellDof<symmTensor>, DgCellDofTenEvalSub,    DgCellDofTenMutSubA, DgCellDofTenMutSubB)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&,  cellDof<tensor>,     cellDof<symmTensor>, cellDof<symmTensor>, DgCellDofTenEvalAnd,    Foam::dgCellDofTmp::NoMutateLeft, Foam::dgCellDofTmp::NoMutateRight)
DG_CELL_DOF_TMP_DEFINE_BINARY_OPERATOR(&&, cellDof<scalar>,     cellDof<symmTensor>, cellDof<symmTensor>, DgCellDofTenEvalAndAnd, Foam::dgCellDofTmp::NoMutateLeft, Foam::dgCellDofTmp::NoMutateRight)

} // End namespace Foam

#endif

// ************************************************************************* //

