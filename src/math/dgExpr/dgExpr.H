/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2026 Ha Nam Tran
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
    Expression-template algebra for GaussField evaluation.

    This header keeps the dgExpr path explicit and opt-in:

        auto expr =
            dg::expr(rhoUG)*dg::expr(UG)
          + dg::expr(pG)*symmTensor::I;

    The important contract is unchanged:

    - GaussField is still the real storage.
    - dgExpr nodes are only lightweight recipes.
    - The recipe is fused only when dg::assign() or dg::evaluate() is called.

    Unlike the original minimal prototype, this version covers the common
    GaussField algebra already exposed by dgGaussFieldMath:

    - binary operators: +, -, *, /, ^, &, &&
    - unary operators/functions: -, mag, magSqr, sqrt, dev, symm, skew,
      T, tr, det, inv
    - literal/uniform mixing: expr(a) + 1.0, 0.5*expr(U), expr(p)*symmTensor::I
    - pow(base, exponent)

    What it still does NOT do on purpose:

    - It does not implicitly convert GaussField to an expression node.
      Users must still write dg::expr(field) to enter the fused path.
    - It does not replace the legacy GaussField operators.
      Both lanes continue to coexist safely.

    Mental model for readers new to expression templates:

    1. A \c GaussField / \c cellGaussField / \c faceGaussField is the real
       storage. It owns the arrays of values on which DG algebra operates.

    2. \c dg::expr(field) does not evaluate anything. It only wraps the field
       in a tiny "leaf" object that knows how to read values from the field.

    3. Operators such as \c +, \c *, \c mag, \c pow do not create a new field.
       They build a tree of small node objects that describe
       "how to compute one value when someone asks for it".

    4. The tree stays lazy until \c dg::assign(dst, expr) or
       \c dg::evaluate(expr, ...) is called. Only at that moment do we loop
       over Gauss points and write concrete values into storage.

    Example:

        auto expr1 = dg::expr(rhoG);
        auto expr2 = dg::expr(UG);
        auto flux  = expr1*expr2 + dg::expr(pG)*symmTensor::I;

    After the three lines above:

    - \c rhoG, \c UG and \c pG still own all numerical data.
    - \c flux owns no field data of its own.
    - \c flux only stores a compact recipe/tree that references those inputs.

    This is the key reason for the subsystem:

    - users can write readable algebra;
    - intermediate fields are avoided;
    - evaluation is fused into the final assignment/evaluate pass.

    Temporary objects follow the same idea, but must be kept alive long enough.
    For that reason the subsystem has dedicated "tmp leaf" wrappers which own
    the \c tmp<T> handle when an expression is built from a temporary field.

\*---------------------------------------------------------------------------*/

/**
 * \file
 * \brief Umbrella include for the DG expression-template subsystem.
 *
 * Include this header when user code needs the whole expression-template API:
 *
 * - leaf adapters (\c expr(field), \c uniform(value))
 * - algebraic node builders (operators/functions)
 * - comparison helpers
 * - evaluation/materialization helpers
 */

#ifndef Foam_dgExpr_H
#define Foam_dgExpr_H

#include "GaussField.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "tmp.H"
#include <type_traits>
#include <utility>

#include "../dgExpr/dgExprTraits.H"
#include "../dgExpr/dgExprLeaves.H"
#include "../dgExpr/dgExprOps.H"
#include "../dgExpr/dgExprNodes.H"
#include "../dgExpr/dgExprFunctions.H"
#include "../dgExpr/dgExprCompare.H"
#include "../dgExpr/dgExprEvaluate.H"

#endif

// ************************************************************************* //
