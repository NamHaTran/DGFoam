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
 *
 * The important teaching point in this file is that the "recipe" of an
 * expression is literally stored in the data members of these classes:
 *
 * - a unary node stores one child object;
 * - a binary node stores two child objects;
 * - the template parameter \c Op says how to combine child values.
 *
 * So for a user expression like
 *
 *     dg::expr(rhoU)/dg::expr(rho)
 *
 * the instantiated type is roughly
 *
 *     binary_expr<div_op, gauss_field_ref<vector>, gauss_field_ref<scalar>>
 *
 * and the recipe lives inside the members \c lhs_ and \c rhs_.
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
 * Template parameters:
 *
 * - \c Op: the operation to apply pointwise, such as \c neg_op or \c mag_op;
 * - \c Expr: the child node type whose values feed that operation.
 *
 * Stores exactly one child expression. The result at quadrature point \c gpI
 * is computed on demand as:
 *
 *     Op::apply(child.valueAt(gpI))
 *
 * No intermediate field is materialized.
 *
 * In other words, this class does not "hold computed values". It holds the
 * instruction "when asked for one point, ask my child for one point and then
 * transform it with \c Op".
 */
template<class Op, class Expr>
class unary_expr
{
    /// \brief Child recipe node from which pointwise input values are read.
    Expr operand_;

public:

    using expression_tag = void;
    static constexpr bool supports_cell = expr_supports_cell_v<Expr>;
    static constexpr bool supports_face = expr_supports_face_v<Expr>;
    static constexpr expr_domain_kind domain = expr_domain_v<Expr>;

    /**
     * \brief Pointwise result type of this node.
     *
     * The type is deduced by pretending to call \c Op::apply() on one sample
     * value produced by the child expression.
     */
    using value_type = std::decay_t
    <
        decltype(Op::apply(std::declval<expr_result_t<Expr>>()))
    >;

    /// \brief Store the child recipe by value.
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
 *
 * The recipe here is:
 *
 * - keep a callable in \c func_;
 * - keep one child node in \c operand_;
 * - when queried, apply \c func_ to one value produced by \c operand_.
 */
template<class Func, class Expr>
class mapped_unary_expr
{
    /// \brief User-provided pointwise transformation.
    Func func_;
    /// \brief Child recipe node.
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
 * Template parameters:
 *
 * - \c Op: the pointwise binary operation such as \c add_op or \c mul_op;
 * - \c LExpr: left child node type;
 * - \c RExpr: right child node type.
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
 *
 * This class is where the expression tree really becomes visible:
 *
 * - \c lhs_ stores the left recipe subtree;
 * - \c rhs_ stores the right recipe subtree;
 * - \c Op determines how the two sampled values are combined.
 */
template<class Op, class LExpr, class RExpr>
class binary_expr
{
    /// \brief Left child recipe subtree.
    LExpr lhs_;
    /// \brief Right child recipe subtree.
    RExpr rhs_;

public:

    using expression_tag = void;
    static constexpr bool supports_cell =
        expr_supports_cell_v<LExpr> && expr_supports_cell_v<RExpr>;
    static constexpr bool supports_face =
        expr_supports_face_v<LExpr> && expr_supports_face_v<RExpr>;
    static constexpr expr_domain_kind domain =
        merge_expr_domain(expr_domain_v<LExpr>, expr_domain_v<RExpr>);

    /**
     * \brief Pointwise result type of this binary node.
     *
     * The type is deduced from one hypothetical application of \c Op to one
     * value from the left child and one value from the right child.
     */
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

    /// \brief Store both child recipe subtrees by value.
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
 *
 * It behaves like \c binary_expr, except that the operation is stored as a
 * runtime callable in \c func_ instead of being encoded purely as a type.
 */
template<class Func, class LExpr, class RExpr>
class mapped_binary_expr
{
    /// \brief User-provided pointwise binary operation.
    Func func_;
    /// \brief Left child recipe subtree.
    LExpr lhs_;
    /// \brief Right child recipe subtree.
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
