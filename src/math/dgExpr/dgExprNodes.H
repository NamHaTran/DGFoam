/**
 * \file
 * \brief Internal node types that build lazy DG expression trees.
 *
 * These classes are the "recipe objects" of the subsystem.
 *
 * They deliberately store child expressions, not evaluated DG fields. For
 * example, the expression
 *
 *     dg::expr(a) + 2.0*dg::expr(b)
 *
 * becomes a small tree of node objects. Each node knows how to compute a
 * single quadrature-point result when asked, but no loop over all points
 * happens here. The actual loop is deferred to \c dg::assign() or
 * \c dg::evaluate().
 */

#ifndef Foam_dgExprNodes_H
#define Foam_dgExprNodes_H

namespace Foam
{
namespace dg
{

/**
 * \brief Unary expression node backed by a stateless operator functor.
 *
 * Stores exactly one child expression. The result at quadrature point \c gpI
 * is computed on demand as:
 *
 *     Op::apply(child.valueAt(gpI))
 *
 * No intermediate field is materialized.
 */
template<class Op, class Expr>
class unary_expr
{
    Expr operand_;

public:

    using expression_tag = void;
    static constexpr bool supports_cell = expr_supports_cell_v<Expr>;
    static constexpr bool supports_face = expr_supports_face_v<Expr>;
    static constexpr expr_domain_kind domain = expr_domain_v<Expr>;

    using value_type = std::decay_t
    <
        decltype(Op::apply(std::declval<expr_result_t<Expr>>()))
    >;

    explicit unary_expr(const Expr& operand)
    :
        operand_(operand)
    {}

    inline label cellSizeHint() const
    {
        return operand_.cellSizeHint();
    }

    inline label faceSizeHint() const
    {
        return operand_.faceSizeHint();
    }

    inline value_type cellValue(const label gpI) const
    {
        return Op::apply(operand_.cellValue(gpI));
    }

    inline value_type faceMinusValue(const label gpI) const
    {
        return Op::apply(operand_.faceMinusValue(gpI));
    }

    inline value_type facePlusValue(const label gpI) const
    {
        return Op::apply(operand_.facePlusValue(gpI));
    }

    inline value_type faceMinusValueOnFace(const label faceI, const label gpI) const
    {
        return Op::apply(operand_.faceMinusValueOnFace(faceI, gpI));
    }

    inline value_type facePlusValueOnFace(const label faceI, const label gpI) const
    {
        return Op::apply(operand_.facePlusValueOnFace(faceI, gpI));
    }
};


/**
 * \brief Unary expression node backed by a user-supplied callable.
 *
 * This is the customizable version of \c unary_expr. It stores the callable by
 * value together with the operand, which makes it suitable for small lambdas
 * and function objects used through \c dg::map().
 */
template<class Func, class Expr>
class mapped_unary_expr
{
    Func func_;
    Expr operand_;

public:

    using expression_tag = void;
    static constexpr bool supports_cell = expr_supports_cell_v<Expr>;
    static constexpr bool supports_face = expr_supports_face_v<Expr>;
    static constexpr expr_domain_kind domain = expr_domain_v<Expr>;

    using value_type = std::decay_t
    <
        decltype
        (
            std::declval<const Func&>()(std::declval<expr_result_t<Expr>>())
        )
    >;

    mapped_unary_expr(const Func& func, const Expr& operand)
    :
        func_(func),
        operand_(operand)
    {}

    inline label cellSizeHint() const
    {
        return operand_.cellSizeHint();
    }

    inline label faceSizeHint() const
    {
        return operand_.faceSizeHint();
    }

    inline value_type cellValue(const label gpI) const
    {
        return func_(operand_.cellValue(gpI));
    }

    inline value_type faceMinusValue(const label gpI) const
    {
        return func_(operand_.faceMinusValue(gpI));
    }

    inline value_type facePlusValue(const label gpI) const
    {
        return func_(operand_.facePlusValue(gpI));
    }

    inline value_type faceMinusValueOnFace(const label faceI, const label gpI) const
    {
        return func_(operand_.faceMinusValueOnFace(faceI, gpI));
    }

    inline value_type facePlusValueOnFace(const label faceI, const label gpI) const
    {
        return func_(operand_.facePlusValueOnFace(faceI, gpI));
    }
};


/**
 * \brief Binary expression node backed by a stateless operator functor.
 *
 * A binary node stores its left and right children by value. Those child
 * objects are themselves usually just tiny references or other nodes, not full
 * DG fields. Reading the tree recursively shows the original algebra.
 *
 * Example:
 *
 *     auto e = dg::expr(rho)*dg::expr(U);
 *
 * creates one \c binary_expr<mul_op, ...> holding two leaves. The product is
 * only computed when \c cellValue() / \c faceValue() is queried during the
 * final evaluation pass.
 */
template<class Op, class LExpr, class RExpr>
class binary_expr
{
    LExpr lhs_;
    RExpr rhs_;

public:

    using expression_tag = void;
    static constexpr bool supports_cell =
        expr_supports_cell_v<LExpr> && expr_supports_cell_v<RExpr>;
    static constexpr bool supports_face =
        expr_supports_face_v<LExpr> && expr_supports_face_v<RExpr>;
    static constexpr expr_domain_kind domain =
        merge_expr_domain(expr_domain_v<LExpr>, expr_domain_v<RExpr>);

    using value_type = std::decay_t
    <
        decltype
        (
            Op::apply
            (
                std::declval<expr_result_t<LExpr>>(),
                std::declval<expr_result_t<RExpr>>()
            )
        )
    >;

    binary_expr(const LExpr& lhs, const RExpr& rhs)
    :
        lhs_(lhs),
        rhs_(rhs)
    {}

    inline label cellSizeHint() const
    {
        return merge_expr_extent(lhs_.cellSizeHint(), rhs_.cellSizeHint());
    }

    inline label faceSizeHint() const
    {
        return merge_expr_extent(lhs_.faceSizeHint(), rhs_.faceSizeHint());
    }

    inline value_type cellValue(const label gpI) const
    {
        return Op::apply(lhs_.cellValue(gpI), rhs_.cellValue(gpI));
    }

    inline value_type faceMinusValue(const label gpI) const
    {
        return Op::apply(lhs_.faceMinusValue(gpI), rhs_.faceMinusValue(gpI));
    }

    inline value_type facePlusValue(const label gpI) const
    {
        return Op::apply(lhs_.facePlusValue(gpI), rhs_.facePlusValue(gpI));
    }

    inline value_type faceMinusValueOnFace(const label faceI, const label gpI) const
    {
        return Op::apply
        (
            lhs_.faceMinusValueOnFace(faceI, gpI),
            rhs_.faceMinusValueOnFace(faceI, gpI)
        );
    }

    inline value_type facePlusValueOnFace(const label faceI, const label gpI) const
    {
        return Op::apply
        (
            lhs_.facePlusValueOnFace(faceI, gpI),
            rhs_.facePlusValueOnFace(faceI, gpI)
        );
    }
};


/**
 * \brief Binary expression node backed by a user-supplied callable.
 *
 * Used by \c dg::map(func, lhs, rhs) for pointwise custom operations while
 * preserving the same lazy evaluation model.
 */
template<class Func, class LExpr, class RExpr>
class mapped_binary_expr
{
    Func func_;
    LExpr lhs_;
    RExpr rhs_;

public:

    using expression_tag = void;
    static constexpr bool supports_cell =
        expr_supports_cell_v<LExpr> && expr_supports_cell_v<RExpr>;
    static constexpr bool supports_face =
        expr_supports_face_v<LExpr> && expr_supports_face_v<RExpr>;
    static constexpr expr_domain_kind domain =
        merge_expr_domain(expr_domain_v<LExpr>, expr_domain_v<RExpr>);

    using value_type = std::decay_t
    <
        decltype
        (
            std::declval<const Func&>()(
                std::declval<expr_result_t<LExpr>>(),
                std::declval<expr_result_t<RExpr>>()
            )
        )
    >;

    mapped_binary_expr(const Func& func, const LExpr& lhs, const RExpr& rhs)
    :
        func_(func),
        lhs_(lhs),
        rhs_(rhs)
    {}

    inline label cellSizeHint() const
    {
        return merge_expr_extent(lhs_.cellSizeHint(), rhs_.cellSizeHint());
    }

    inline label faceSizeHint() const
    {
        return merge_expr_extent(lhs_.faceSizeHint(), rhs_.faceSizeHint());
    }

    inline value_type cellValue(const label gpI) const
    {
        return func_(lhs_.cellValue(gpI), rhs_.cellValue(gpI));
    }

    inline value_type faceMinusValue(const label gpI) const
    {
        return func_(lhs_.faceMinusValue(gpI), rhs_.faceMinusValue(gpI));
    }

    inline value_type facePlusValue(const label gpI) const
    {
        return func_(lhs_.facePlusValue(gpI), rhs_.facePlusValue(gpI));
    }

    inline value_type faceMinusValueOnFace(const label faceI, const label gpI) const
    {
        return func_
        (
            lhs_.faceMinusValueOnFace(faceI, gpI),
            rhs_.faceMinusValueOnFace(faceI, gpI)
        );
    }

    inline value_type facePlusValueOnFace(const label faceI, const label gpI) const
    {
        return func_
        (
            lhs_.facePlusValueOnFace(faceI, gpI),
            rhs_.facePlusValueOnFace(faceI, gpI)
        );
    }
};

} // End namespace dg
} // End namespace Foam

#endif
