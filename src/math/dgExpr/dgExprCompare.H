/**
 * \file
 * \brief Pointwise comparison support for DG expression trees.
 *
 * This header is intentionally different in spirit from the algebraic builders.
 *
 * Operators such as \c + or \c * build a new lazy recipe node.
 * Comparison operators such as \c < or \c == do not build another DG
 * expression. They immediately reduce the pointwise comparison over all
 * supported quadrature points and return one \c bool.
 *
 * So:
 *
 *     dg::expr(a) + dg::expr(b)
 *
 * stays lazy, whereas
 *
 *     dg::expr(a) < dg::expr(b)
 *
 * is eager and produces a truth value right away.
 */

#ifndef Foam_dgExprCompare_H
#define Foam_dgExprCompare_H

namespace Foam
{
namespace dg
{

/**
 * \brief Detect whether a predicate may compare the results of two operands.
 *
 * This trait checks two things:
 *
 * - whether \c Pred::apply(lhsValue, rhsValue) is a valid expression;
 * - whether its result is convertible to \c bool.
 *
 * It keeps the comparison API available only for sensible value types.
 */
template<class Pred, class LExpr, class RExpr, class = void>
struct is_expr_comparison
:
    std::false_type
{};


template<class Pred, class LExpr, class RExpr>
struct is_expr_comparison
<
    Pred,
    LExpr,
    RExpr,
    std::void_t
    <
        decltype
        (
            Pred::apply
            (
                std::declval<expr_result_t<expression_node_t<LExpr>>>(),
                std::declval<expr_result_t<expression_node_t<RExpr>>>()
            )
        )
    >
>
:
    std::is_convertible
    <
        decltype
        (
            Pred::apply
            (
                std::declval<expr_result_t<expression_node_t<LExpr>>>(),
                std::declval<expr_result_t<expression_node_t<RExpr>>>()
            )
        ),
        bool
    >
{};


template<class Pred, class LExpr, class RExpr>
inline constexpr bool is_expr_comparison_v =
    is_expr_comparison<Pred, LExpr, RExpr>::value;

/// Pointwise \c < predicate.
struct less_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs < rhs;
    }
};


/// Pointwise \c <= predicate.
struct less_equal_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs <= rhs;
    }
};


/// Pointwise \c > predicate.
struct greater_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs > rhs;
    }
};


/// Pointwise \c >= predicate.
struct greater_equal_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs >= rhs;
    }
};


/// Pointwise \c == predicate.
struct equal_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs == rhs;
    }
};


/// Pointwise \c != predicate.
struct not_equal_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs != rhs;
    }
};


/**
 * \brief Compare two normalized expression nodes over all supported values.
 *
 * The comparison is reduced with logical-and semantics:
 *
 * - every sampled point must satisfy the predicate;
 * - the first failing point returns \c false immediately;
 * - if every reachable point passes, the final answer is \c true.
 *
 * This makes the operators read naturally as "all pointwise values satisfy the
 * relation".
 */
template<class Pred, class LNode, class RNode>
inline bool compare_nodes(const LNode& lhs, const RNode& rhs)
{
    bool compared = false;

    if constexpr (expr_supports_cell_v<LNode> && expr_supports_cell_v<RNode>)
    {
        const label nCell = merge_expr_extent(lhs.cellSizeHint(), rhs.cellSizeHint());

        if (nCell != unknown_expr_extent)
        {
            compared = true;

            for (label gpI = 0; gpI < nCell; ++gpI)
            {
                if (!Pred::apply(lhs.cellValue(gpI), rhs.cellValue(gpI)))
                {
                    return false;
                }
            }
        }
    }

    if constexpr (expr_supports_face_v<LNode> && expr_supports_face_v<RNode>)
    {
        const label nFace = merge_expr_extent(lhs.faceSizeHint(), rhs.faceSizeHint());

        if (nFace != unknown_expr_extent)
        {
            compared = true;

            for (label gpI = 0; gpI < nFace; ++gpI)
            {
                if (!Pred::apply(lhs.faceMinusValue(gpI), rhs.faceMinusValue(gpI)))
                {
                    return false;
                }

                if (!Pred::apply(lhs.facePlusValue(gpI), rhs.facePlusValue(gpI)))
                {
                    return false;
                }
            }
        }
    }

    if (!compared)
    {
        if constexpr (expr_supports_cell_v<LNode> && expr_supports_cell_v<RNode>)
        {
            return Pred::apply(lhs.cellValue(0), rhs.cellValue(0));
        }
        else if constexpr (expr_supports_face_v<LNode> && expr_supports_face_v<RNode>)
        {
            return
                Pred::apply(lhs.faceMinusValue(0), rhs.faceMinusValue(0))
             && Pred::apply(lhs.facePlusValue(0), rhs.facePlusValue(0));
        }
    }

    return true;
}


/**
 * \brief Compare two expression operands after normalizing literals to
 *        uniform-expression leaves.
 *
 * This is the comparison counterpart of the algebraic builders in
 * \c dgExprFunctions.H: user operands are normalized first, then handed to the
 * internal implementation.
 */
template<class Pred, class LExpr, class RExpr>
inline auto compare_expressions(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>
     && is_expr_comparison_v<Pred, LExpr, RExpr>,
        bool
    >
{
    return compare_nodes<Pred>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise less-than reduction across two DG expression operands.
template<class LExpr, class RExpr>
inline auto operator<(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>
     && is_expr_comparison_v<less_op, LExpr, RExpr>,
        bool
    >
{
    return compare_expressions<less_op>(lhs, rhs);
}


/// Pointwise less-than-or-equal reduction across two DG expression operands.
template<class LExpr, class RExpr>
inline auto operator<=(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>
     && is_expr_comparison_v<less_equal_op, LExpr, RExpr>,
        bool
    >
{
    return compare_expressions<less_equal_op>(lhs, rhs);
}


/// Pointwise greater-than reduction across two DG expression operands.
template<class LExpr, class RExpr>
inline auto operator>(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>
     && is_expr_comparison_v<greater_op, LExpr, RExpr>,
        bool
    >
{
    return compare_expressions<greater_op>(lhs, rhs);
}


/// Pointwise greater-than-or-equal reduction across two DG expression operands.
template<class LExpr, class RExpr>
inline auto operator>=(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>
     && is_expr_comparison_v<greater_equal_op, LExpr, RExpr>,
        bool
    >
{
    return compare_expressions<greater_equal_op>(lhs, rhs);
}


/// Pointwise equality reduction across two DG expression operands.
template<class LExpr, class RExpr>
inline auto operator==(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>
     && is_expr_comparison_v<equal_op, LExpr, RExpr>,
        bool
    >
{
    return compare_expressions<equal_op>(lhs, rhs);
}


/// Pointwise inequality reduction across two DG expression operands.
template<class LExpr, class RExpr>
inline auto operator!=(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>
     && is_expr_comparison_v<not_equal_op, LExpr, RExpr>,
        bool
    >
{
    return compare_expressions<not_equal_op>(lhs, rhs);
}

} // End namespace dg
} // End namespace Foam

#endif
