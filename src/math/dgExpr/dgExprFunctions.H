/**
 * \file
 * \brief User-facing builders and algebraic helpers for DG expressions.
 *
 * This is the main "how do I create an expression?" entry point.
 *
 * Common user flow:
 *
 *     auto e =
 *         dg::expr(rhoG)*dg::expr(UG)
 *       + dg::expr(pG)*symmTensor::I;
 *
 *     dg::assign(result, e);        // write into existing storage
 *     auto tResult = dg::evaluate(e, cellID, mesh);  // allocate new storage
 *
 * Reading that flow from top to bottom:
 *
 * - \c expr(field) wraps storage into a leaf;
 * - operators/functions build a lazy tree;
 * - \c assign/evaluate are the places where values are finally computed.
 */

#ifndef Foam_dgExprFunctions_H
#define Foam_dgExprFunctions_H

namespace Foam
{
namespace dg
{

/**
 * \brief Enter the lazy DG-expression path from a full \c GaussField.
 *
 * Important: this returns a lightweight leaf object, not a new field.
 */
template<class Type>
inline gauss_field_ref<Type> expr(const GaussField<Type>& field)
{
    return gauss_field_ref<Type>(field);
}


/**
 * \brief Enter the lazy DG-expression path from a \c cellGaussField.
 *
 * Use this when the source already lives only on cell quadrature points.
 */
template<class Type>
inline cell_gauss_field_ref<Type> expr(const cellGaussField<Type>& field)
{
    return cell_gauss_field_ref<Type>(field);
}


/**
 * \brief Enter the lazy DG-expression path from a \c faceGaussField.
 *
 * Use this when the source already lives only on face quadrature points.
 */
template<class Type>
inline face_gauss_field_ref<Type> expr(const faceGaussField<Type>& field)
{
    return face_gauss_field_ref<Type>(field);
}


/**
 * \brief Wrap a temporary full \c GaussField as an owning expression leaf.
 *
 * This preserves the lifetime of temporary storage until the whole expression
 * tree is evaluated.
 */
template<class Type>
inline tmp_gauss_field_ref<Type> expr(const tmp<GaussField<Type>>& field)
{
    return tmp_gauss_field_ref<Type>(field);
}


/**
 * \brief Wrap a temporary \c cellGaussField as an owning expression leaf.
 *
 * Same idea as the full-field overload, but restricted to cell storage.
 */
template<class Type>
inline tmp_cell_gauss_field_ref<Type> expr(const tmp<cellGaussField<Type>>& field)
{
    return tmp_cell_gauss_field_ref<Type>(field);
}


/**
 * \brief Wrap a temporary \c faceGaussField as an owning expression leaf.
 *
 * Same idea as the full-field overload, but restricted to face storage.
 */
template<class Type>
inline tmp_face_gauss_field_ref<Type> expr(const tmp<faceGaussField<Type>>& field)
{
    return tmp_face_gauss_field_ref<Type>(field);
}


/**
 * \brief Promote a literal value to a uniform DG expression.
 *
 * This is the storage model for constants inside an expression tree. The value
 * is copied once into the leaf and then read identically at every quadrature
 * point.
 */
template<class Type>
inline auto uniform(const Type& value)
    -> std::enable_if_t<is_uniform_value_v<Type>, uniform_expr<std::decay_t<Type>>>
{
    return uniform_expr<std::decay_t<Type>>(value);
}


/**
 * \brief Apply a callable pointwise to a single expression.
 *
 * Example:
 *
 *     auto positive = dg::map([](const scalar x) { return max(x, 0.0); }, e);
 */
template<class Func, class Expr>
inline auto map(const Func& func, const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        mapped_unary_expr<std::decay_t<Func>, std::decay_t<Expr>>
    >
{
    return mapped_unary_expr<std::decay_t<Func>, std::decay_t<Expr>>
    (
        func,
        expression
    );
}


/**
 * \brief Apply a callable pointwise to two expressions or expression operands.
 *
 * Literals are accepted as long as at least one side is already an expression.
 * Non-expression operands are normalized into uniform leaves internally.
 */
template<class Func, class LExpr, class RExpr>
inline auto map(const Func& func, const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        mapped_binary_expr
        <
            std::decay_t<Func>,
            expression_node_t<LExpr>,
            expression_node_t<RExpr>
        >
    >
{
    return mapped_binary_expr
    <
        std::decay_t<Func>,
        expression_node_t<LExpr>,
        expression_node_t<RExpr>
    >
    (
        func,
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Unary minus on an expression tree. Returns a new lazy node, not a field.
template<class Expr>
inline auto operator-(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<neg_op, std::decay_t<Expr>>
    >
{
    return unary_expr<neg_op, std::decay_t<Expr>>(expression);
}


/// Pointwise addition between DG expression operands. Still lazy.
template<class LExpr, class RExpr>
inline auto operator+(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr<add_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    >
{
    return binary_expr<add_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise subtraction between DG expression operands. Still lazy.
template<class LExpr, class RExpr>
inline auto operator-(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr<sub_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    >
{
    return binary_expr<sub_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise multiplication between DG expression operands. Still lazy.
template<class LExpr, class RExpr>
inline auto operator*(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr<mul_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    >
{
    return binary_expr<mul_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise division between DG expression operands. Still lazy.
template<class LExpr, class RExpr>
inline auto operator/(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr<div_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    >
{
    return binary_expr<div_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise \c operator^ between DG expression operands.
template<class LExpr, class RExpr>
inline auto operator^(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr<cross_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    >
{
    return binary_expr<cross_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise \c operator& between DG expression operands.
template<class LExpr, class RExpr>
inline auto operator&(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr<dot_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    >
{
    return binary_expr<dot_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise \c operator&& between DG expression operands.
template<class LExpr, class RExpr>
inline auto operator&&(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr
        <
            double_dot_op,
            expression_node_t<LExpr>,
            expression_node_t<RExpr>
        >
    >
{
    return binary_expr
    <
        double_dot_op,
        expression_node_t<LExpr>,
        expression_node_t<RExpr>
    >
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/**
 * \brief Pointwise power between DG expression operands.
 */
template<class LExpr, class RExpr>
inline auto pow(const LExpr& lhs, const RExpr& rhs)
    -> std::enable_if_t
    <
        (is_expression_v<LExpr> || is_expression_v<RExpr>)
     && is_expr_operand_v<LExpr>
     && is_expr_operand_v<RExpr>
     && is_expr_domain_compatible_v<expression_node_t<LExpr>, expression_node_t<RExpr>>,
        binary_expr<pow_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    >
{
    return binary_expr<pow_op, expression_node_t<LExpr>, expression_node_t<RExpr>>
    (
        make_expression_node(lhs),
        make_expression_node(rhs)
    );
}


/// Pointwise magnitude.
template<class Expr>
inline auto mag(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<mag_op, std::decay_t<Expr>>
    >
{
    return unary_expr<mag_op, std::decay_t<Expr>>(expression);
}


/// Pointwise squared magnitude.
template<class Expr>
inline auto magSqr(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<mag_sqr_op, std::decay_t<Expr>>
    >
{
    return unary_expr<mag_sqr_op, std::decay_t<Expr>>(expression);
}


/// Pointwise square root.
template<class Expr>
inline auto sqrt(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<sqrt_op, std::decay_t<Expr>>
    >
{
    return unary_expr<sqrt_op, std::decay_t<Expr>>(expression);
}


/// Pointwise deviatoric part.
template<class Expr>
inline auto dev(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<dev_op, std::decay_t<Expr>>
    >
{
    return unary_expr<dev_op, std::decay_t<Expr>>(expression);
}


/// Pointwise symmetric part.
template<class Expr>
inline auto symm(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<symm_op, std::decay_t<Expr>>
    >
{
    return unary_expr<symm_op, std::decay_t<Expr>>(expression);
}


/// Pointwise skew-symmetric part.
template<class Expr>
inline auto skew(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<skew_op, std::decay_t<Expr>>
    >
{
    return unary_expr<skew_op, std::decay_t<Expr>>(expression);
}


/// Pointwise transpose.
template<class Expr>
inline auto T(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<transpose_op, std::decay_t<Expr>>
    >
{
    return unary_expr<transpose_op, std::decay_t<Expr>>(expression);
}


/// Pointwise trace.
template<class Expr>
inline auto tr(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<trace_op, std::decay_t<Expr>>
    >
{
    return unary_expr<trace_op, std::decay_t<Expr>>(expression);
}


/// Pointwise determinant.
template<class Expr>
inline auto det(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<det_op, std::decay_t<Expr>>
    >
{
    return unary_expr<det_op, std::decay_t<Expr>>(expression);
}


/// Pointwise inverse.
template<class Expr>
inline auto inv(const Expr& expression)
    -> std::enable_if_t
    <
        is_expression_v<Expr>,
        unary_expr<inv_op, std::decay_t<Expr>>
    >
{
    return unary_expr<inv_op, std::decay_t<Expr>>(expression);
}

} // End namespace dg
} // End namespace Foam

#endif
