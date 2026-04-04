/**
 * \file
 * \brief Evaluation, reduction and assignment helpers for DG expressions.
 *
 * This file is the boundary between "lazy recipe" and "real DG storage".
 *
 * Everything before this point builds descriptions of computations. The
 * functions below are where we finally walk the expression tree over all
 * quadrature points and either:
 *
 * - reduce the values to a scalar/aggregate result, or
 * - write them into concrete DG fields.
 *
 * This is the execution boundary of \c dgExpr. If earlier headers answer
 * "what is the recipe?", this header answers "when and how is that recipe run?".
 */

#ifndef Foam_dgExprEvaluate_H
#define Foam_dgExprEvaluate_H

namespace Foam
{
namespace dg
{

namespace detail
{

/**
 * \brief Visit every runtime value reachable from an expression.
 *
 * \return True when the expression provided a known runtime extent and was
 *         fully reduced through iteration.
 *
 * This helper is used by reductions such as \c maxValue, \c minValue, \c sum
 * and \c average. It centralizes the "walk every accessible DG point" logic so
 * those algorithms only need to describe what to do with each sampled value.
 *
 * Conceptually, this is the reduction analogue of \c assign():
 *
 * - ask the expression whether it supports cell and/or face access;
 * - ask for the number of available points;
 * - visit each point exactly once;
 * - let the caller decide what to do with each sampled value.
 */
template<class Expr, class Func>
inline bool forEachExprValue(const Expr& expression, Func&& func)
{
    bool reduced = false;

    if constexpr (expr_supports_cell_v<Expr>)
    {
        const label nCell = expression.cellSizeHint();

        if (nCell != unknown_expr_extent)
        {
            reduced = true;

            for (label gpI = 0; gpI < nCell; ++gpI)
            {
                func(expression.cellValue(gpI));
            }
        }
    }

    if constexpr (expr_supports_face_v<Expr>)
    {
        const label nFace = expression.faceSizeHint();

        if (nFace != unknown_expr_extent)
        {
            reduced = true;

            for (label gpI = 0; gpI < nFace; ++gpI)
            {
                func(expression.faceMinusValue(gpI));
                func(expression.facePlusValue(gpI));
            }
        }
    }

    return reduced;
}


/**
 * \brief Return a representative value for expressions with no known extent.
 *
 * Uniform expressions may have no natural runtime size. In that case we still
 * need one valid sample value for reductions and comparisons, and any
 * quadrature point gives the same answer.
 */
template<class Expr>
inline expr_result_t<Expr> firstExprValue(const Expr& expression)
{
    if constexpr (expr_supports_cell_v<Expr>)
    {
        return expression.cellValue(0);
    }
    else
    {
        return expression.faceMinusValue(0);
    }
}

} // End namespace detail


/**
 * \brief Compute the maximum scalar value stored in an expression.
 *
 * This is an eager reduction. Unlike the algebraic builders, it returns an
 * actual scalar immediately.
 */
template<class Expr>
inline std::enable_if_t
<
    is_expression_v<Expr>
 && std::is_same<std::decay_t<expr_result_t<Expr>>, scalar>::value,
    scalar
>
maxValue(const Expr& expression)
{
    scalar result = -GREAT;
    const bool reduced = detail::forEachExprValue
    (
        expression,
        [&](const scalar value)
        {
            result = max(result, value);
        }
    );

    if (!reduced)
    {
        return detail::firstExprValue(expression);
    }

    return result;
}


/**
 * \brief Compute the minimum scalar value stored in an expression.
 */
template<class Expr>
inline std::enable_if_t
<
    is_expression_v<Expr>
 && std::is_same<std::decay_t<expr_result_t<Expr>>, scalar>::value,
    scalar
>
minValue(const Expr& expression)
{
    scalar result = GREAT;
    const bool reduced = detail::forEachExprValue
    (
        expression,
        [&](const scalar value)
        {
            result = min(result, value);
        }
    );

    if (!reduced)
    {
        return detail::firstExprValue(expression);
    }

    return result;
}


/**
 * \brief Sum all values carried by an expression.
 */
template<class Expr>
inline std::enable_if_t<is_expression_v<Expr>, expr_result_t<Expr>>
sum(const Expr& expression)
{
    using value_type = expr_result_t<Expr>;

    value_type result(pTraits<value_type>::zero);

    const bool reduced = detail::forEachExprValue
    (
        expression,
        [&](const value_type& value)
        {
            result += value;
        }
    );

    if (!reduced)
    {
        return detail::firstExprValue(expression);
    }

    return result;
}


/**
 * \brief Compute the average of all values carried by an expression.
 */
template<class Expr>
inline auto average(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        std::decay_t<decltype(std::declval<expr_result_t<Expr>>() / scalar(1))>
    >
{
    using value_type = expr_result_t<Expr>;
    value_type result(pTraits<value_type>::zero);
    label count = 0;

    const bool reduced = detail::forEachExprValue
    (
        expression,
        [&](const value_type& value)
        {
            result += value;
            ++count;
        }
    );

    if (count > 0)
    {
        return result/scalar(count);
    }

    if (!reduced)
    {
        return detail::firstExprValue(expression);
    }

    WarningInFunction
        << "empty expression, returning zero" << endl;

    return pTraits<value_type>::zero;
}


/**
 * \brief Evaluate an expression directly into a \c cellGaussField.
 *
 * This is the fused execution step for cell-only storage. The destination
 * field already owns memory; we simply loop once and ask the expression tree
 * for the value at each quadrature point.
 *
 * Use \c assign(dst, expr) when:
 *
 * - storage already exists and should be reused;
 * - the caller wants to avoid a fresh allocation;
 * - the natural execution model is "write into this field".
 */
template<class Type, class Expr>
inline std::enable_if_t<is_expression_v<Expr>, void> assign
(
    cellGaussField<Type>& dst,
    const Expr& expression
)
{
    static_assert
    (
        std::is_convertible<expr_result_t<Expr>, Type>::value,
        "dg::assign(dst, expr): destination cellGaussField type is "
        "incompatible with the expression result type"
    );

    for (label gpI = 0; gpI < dst.size(); ++gpI)
    {
        dst[gpI] = expression.cellValue(gpI);
    }
}


/**
 * \brief Evaluate an expression directly into a \c faceGaussField.
 *
 * The lazy tree is traversed pointwise and written into the existing face
 * storage, separately for minus and plus traces.
 */
template<class Type, class Expr>
inline std::enable_if_t<is_expression_v<Expr>, void> assign
(
    faceGaussField<Type>& dst,
    const Expr& expression
)
{
    static_assert
    (
        std::is_convertible<expr_result_t<Expr>, Type>::value,
        "dg::assign(dst, expr): destination faceGaussField type is "
        "incompatible with the expression result type"
    );

    for (label gpI = 0; gpI < dst.nGauss(); ++gpI)
    {
        dst.minusValueAt(gpI) = expression.faceMinusValue(gpI);
        dst.plusValueAt(gpI) = expression.facePlusValue(gpI);
    }
}


/**
 * \brief Evaluate an expression directly into a full \c GaussField.
 *
 * This is the most complete materialization path: one fused pass for the cell
 * part and one fused pass for the face part of the destination field.
 */
template<class Type, class Expr>
inline std::enable_if_t<is_expression_v<Expr>, void> assign
(
    GaussField<Type>& dst,
    const Expr& expression
)
{
    static_assert
    (
        std::is_convertible<expr_result_t<Expr>, Type>::value,
        "dg::assign(dst, expr): destination GaussField type is incompatible "
        "with the expression result type"
    );

    for (label gpI = 0; gpI < dst.cellField().size(); ++gpI)
    {
        dst.cellField()[gpI] = expression.cellValue(gpI);
    }

    for (label gpI = 0; gpI < dst.faceField().nGauss(); ++gpI)
    {
        dst.faceField().minusValueAt(gpI) = expression.faceMinusValue(gpI);
        dst.faceField().plusValueAt(gpI) = expression.facePlusValue(gpI);
    }
}


/**
 * \brief Materialize a cell-only expression into a new \c cellGaussField.
 *
 * Use this when you want a fresh owning field instead of assigning into an
 * existing destination.
 *
 * This is the allocation-producing counterpart of \c assign().
 */
template<class ResultType, class Expr>
inline std::enable_if_t<is_expression_v<Expr>, tmp<cellGaussField<ResultType>>>
evaluateCell
(
    const Expr& expression,
    const label cellID,
    const dgGeomMesh* mesh
)
{
    tmp<cellGaussField<ResultType>> tResult
    (
        new cellGaussField<ResultType>(cellID, mesh)
    );
    assign(tResult.ref(), expression);
    return tResult;
}


/**
 * \brief Materialize a cell-only expression using its natural result type.
 *
 * The result type comes from \c expr_result_t<Expr>.
 */
template<class Expr>
inline std::enable_if_t<is_expression_v<Expr>, tmp<cellGaussField<expr_result_t<Expr>>>>
evaluateCell
(
    const Expr& expression,
    const label cellID,
    const dgGeomMesh* mesh
)
{
    return evaluateCell<expr_result_t<Expr>>(expression, cellID, mesh);
}


/**
 * \brief Materialize a face-only expression into a new \c faceGaussField.
 *
 * This allocates storage first and then delegates to \c assign().
 */
template<class ResultType, class Expr>
inline std::enable_if_t<is_expression_v<Expr>, tmp<faceGaussField<ResultType>>>
evaluateFace
(
    const Expr& expression,
    const label cellID,
    const dgGeomMesh* mesh
)
{
    tmp<faceGaussField<ResultType>> tResult
    (
        new faceGaussField<ResultType>(cellID, mesh)
    );
    assign(tResult.ref(), expression);
    return tResult;
}


/**
 * \brief Materialize a face-only expression using its natural result type.
 *
 * The result type comes from \c expr_result_t<Expr>.
 */
template<class Expr>
inline std::enable_if_t<is_expression_v<Expr>, tmp<faceGaussField<expr_result_t<Expr>>>>
evaluateFace
(
    const Expr& expression,
    const label cellID,
    const dgGeomMesh* mesh
)
{
    return evaluateFace<expr_result_t<Expr>>(expression, cellID, mesh);
}


/**
 * \brief Materialize an expression into a new full \c GaussField.
 *
 * Conceptually this is:
 *
 * 1. allocate destination storage;
 * 2. traverse the lazy expression once;
 * 3. store computed values into the new field.
 */
template<class ResultType, class Expr>
inline std::enable_if_t<is_expression_v<Expr>, tmp<GaussField<ResultType>>> evaluate
(
    const Expr& expression,
    const label cellID,
    const dgGeomMesh* mesh
)
{
    tmp<GaussField<ResultType>> tResult = GaussField<ResultType>::New(cellID, mesh);
    assign(tResult.ref(), expression);
    return tResult;
}


/**
 * \brief Materialize an expression using its natural result type.
 *
 * This is the most ergonomic full-field execution entry point when the caller
 * does not need to force a different destination type.
 */
template<class Expr>
inline std::enable_if_t<is_expression_v<Expr>, tmp<GaussField<expr_result_t<Expr>>>> evaluate
(
    const Expr& expression,
    const label cellID,
    const dgGeomMesh* mesh
)
{
    return evaluate<expr_result_t<Expr>>(expression, cellID, mesh);
}

} // End namespace dg
} // End namespace Foam

#endif
