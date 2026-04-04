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
 *
 * This header is the best place to study the user-facing surface of the
 * subsystem. Most APIs below do one of only three things:
 *
 * - enter the expression world from a DG storage object;
 * - normalize literals so they can participate in a tree;
 * - build a new recipe node representing one algebraic step.
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
 *
 * The returned leaf simply remembers how to read one interior value and one
 * face trace from the original field. This is the explicit opt-in boundary
 * between "plain DG storage" and "fused expression-template algebra".
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
 *
 * Most user code does not need to call this explicitly because mixed operators
 * already normalize literals automatically. It is still useful when the code
 * wants to be very explicit about intent, for example with \c symmTensor::I.
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
 *
 * The callable is stored inside the resulting recipe node and is invoked later
 * during evaluation, once per quadrature point.
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
 *
 * The "at least one side is already an expression" rule prevents this API from
 * becoming a generic replacement for ordinary scalar algebra outside DG code.
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


/**
 * \brief Unary minus on an expression tree.
 *
 * Returns a new lazy node, not a field. The resulting object only says
 * "negate my child value when sampled".
 */
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


/**
 * \brief Pointwise addition between DG expression operands.
 *
 * Still lazy: the return value is a recipe node storing the left and right
 * operands, not a materialized DG field.
 */
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


/// \brief Pointwise subtraction between DG expression operands. Still lazy.
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


/// \brief Pointwise multiplication between DG expression operands. Still lazy.
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


/// \brief Pointwise division between DG expression operands. Still lazy.
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
 *
 * This forwards to \c Foam::pow pointwise and returns another lazy binary
 * recipe node.
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


/// \brief Pointwise magnitude. Builds a lazy unary recipe node.
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


/// \brief Pointwise squared magnitude. Builds a lazy unary recipe node.
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


/// \brief Pointwise square root. Builds a lazy unary recipe node.
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


/// \brief Pointwise deviatoric part. Builds a lazy unary recipe node.
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


/// \brief Pointwise symmetric part. Builds a lazy unary recipe node.
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


/// \brief Pointwise skew-symmetric part. Builds a lazy unary recipe node.
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


/// \brief Pointwise transpose. Builds a lazy unary recipe node.
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


/// \brief Pointwise trace. Builds a lazy unary recipe node.
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


/// \brief Pointwise determinant. Builds a lazy unary recipe node.
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


/// \brief Pointwise inverse. Builds a lazy unary recipe node.
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
