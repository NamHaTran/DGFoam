/**
 * \file
 * \brief Core traits and domain metadata for the \c Foam::dg expression system.
 *
 * This file is the compile-time "policy layer" of \c dgExpr.
 *
 * Nothing here evaluates Gauss-point values directly. Instead, these helpers
 * answer structural questions about types so the rest of the subsystem can
 * decide, at compile time:
 *
 * - whether a type is already an expression node;
 * - whether a literal may be promoted into the expression world;
 * - what value type an expression produces;
 * - whether two operands may appear in the same tree;
 * - whether a node can supply cell values, face values, or both.
 *
 * A practical reading tip:
 *
 * 1. \c is_expression_v and \c expr_result_t identify the "shape" of a node;
 * 2. \c is_uniform_value_v and \c is_expr_operand_v control API admission;
 * 3. \c expr_domain_kind and the merge helpers prevent invalid mixes such as
 *    cell-only expressions with face-only expressions.
 */

#ifndef Foam_dgExprTraits_H
#define Foam_dgExprTraits_H

namespace Foam
{
namespace dg
{

/**
 * \brief Detect whether a type participates in the DG expression system.
 *
 * The convention is intentionally lightweight: a type is considered an
 * expression if it exposes the marker alias
 *
 *     using expression_tag = void;
 *
 * Leaf nodes and internal recipe nodes define that alias, so the subsystem can
 * recognize them without any common base class or virtual dispatch.
 */
template<class T, class = void>
struct is_expression
:
    std::false_type
{};


template<class T>
struct is_expression<T, std::void_t<typename std::decay_t<T>::expression_tag>>
:
    std::true_type
{};


template<class T>
inline constexpr bool is_expression_v = is_expression<T>::value;


/**
 * \brief Value type produced by an expression node.
 *
 * Example:
 *
 * - \c gauss_field_ref<vector> produces \c vector
 * - \c unary_expr<mag_op, gauss_field_ref<vector>> produces \c scalar
 * - \c binary_expr<mul_op, gauss_field_ref<scalar>, gauss_field_ref<vector>>
 *   produces \c vector
 *
 * This alias is used throughout the subsystem to deduce destination types,
 * reduction accumulator types and operator result types.
 */
template<class T>
using expr_result_t = typename std::decay_t<T>::value_type;


/**
 * \brief Detect values that may be promoted to a uniform expression leaf.
 *
 * These are the scalar/tensor-like literals that can appear directly in user
 * algebra, for example:
 *
 *     dg::expr(p) + 1.0
 *     0.5*dg::expr(U)
 *     dg::expr(p)*symmTensor::I
 *
 * Such values are not expression nodes by themselves, but they can be wrapped
 * into \c uniform_expr so the tree sees a uniform value at every Gauss point.
 */
template<class T>
struct is_uniform_value
:
    std::integral_constant
    <
        bool,
        std::is_arithmetic<std::decay_t<T>>::value
     || std::is_same<std::decay_t<T>, vector>::value
     || std::is_same<std::decay_t<T>, tensor>::value
     || std::is_same<std::decay_t<T>, symmTensor>::value
     || std::is_same<std::decay_t<T>, sphericalTensor>::value
    >
{};


template<class T>
inline constexpr bool is_uniform_value_v = is_uniform_value<T>::value;


/**
 * \brief Detect valid operands for expression construction.
 *
 * An operand is accepted by the user-facing builders if it is either:
 *
 * - already an expression node, or
 * - a literal/uniform value that may be normalized into \c uniform_expr.
 *
 * This is the gate used by most operator overloads in \c dgExprFunctions.H.
 */
template<class T>
struct is_expr_operand
:
    std::integral_constant<bool, is_expression_v<T> || is_uniform_value_v<T>>
{};


template<class T>
inline constexpr bool is_expr_operand_v = is_expr_operand<T>::value;


/**
 * \brief Sentinel used when an expression leaf cannot report its extent directly.
 *
 * Uniform expressions are the main example: a literal such as \c 1.0 has no
 * intrinsic runtime size, so its leaf reports "unknown extent". When that
 * leaf is combined with a real field-backed expression, the field-backed side
 * provides the concrete number of quadrature points.
 */
static constexpr label unknown_expr_extent = -1;


/**
 * \brief Logical domain covered by an expression node.
 *
 * This describes which part of DG storage a node can legally serve:
 *
 * - \c agnostic: a domain-free node such as \c uniform_expr
 * - \c cell: cell-interior Gauss values only
 * - \c face: face Gauss values only
 * - \c gauss: full \c GaussField, meaning both cell and face values
 * - \c incompatible: internal marker used when two domains cannot be merged
 *
 * The enum lets the builders reject nonsensical expressions early, instead of
 * discovering the mismatch deep inside evaluation loops.
 */
enum class expr_domain_kind
{
    agnostic,
    cell,
    face,
    gauss,
    incompatible
};


/**
 * \brief Merge two runtime extents and abort when both are known but differ.
 *
 * Typical cases:
 *
 * - field + literal:
 *   the literal reports \c unknown_expr_extent, so the field extent wins;
 * - field + field:
 *   both extents must match;
 * - literal + literal:
 *   the merged result stays \c unknown_expr_extent.
 *
 * This function enforces the assumption that pointwise operations combine
 * values from matching quadrature layouts.
 */
inline label merge_expr_extent(const label lhsExtent, const label rhsExtent)
{
    if
    (
        lhsExtent != unknown_expr_extent
     && rhsExtent != unknown_expr_extent
     && lhsExtent != rhsExtent
    )
    {
        FatalErrorInFunction
            << "Mismatched expression extents: " << lhsExtent
            << " vs " << rhsExtent
            << abort(FatalError);
    }

    return
        lhsExtent != unknown_expr_extent
      ? lhsExtent
      : rhsExtent;
}


/**
 * \brief Convenience accessor for whether an expression can evaluate cell values.
 *
 * The value comes from the node type itself. Because it is compile-time
 * information, later code can use \c if constexpr and completely remove
 * branches that do not apply to a given expression.
 */
template<class Expr>
inline constexpr bool expr_supports_cell_v = std::decay_t<Expr>::supports_cell;


/**
 * \brief Convenience accessor for whether an expression can evaluate face values.
 *
 * Used in the same way as \c expr_supports_cell_v, but for face traces.
 */
template<class Expr>
inline constexpr bool expr_supports_face_v = std::decay_t<Expr>::supports_face;


/**
 * \brief Merge two expression domains into their compatible common domain.
 *
 * The merge rules are:
 *
 * - agnostic + X -> X
 * - X + X -> X
 * - cell + face -> incompatible
 * - gauss + cell -> incompatible
 * - gauss + face -> incompatible
 *
 * In other words, only literals are domain-neutral. Real storage-backed nodes
 * must already agree on which DG region they represent.
 */
constexpr expr_domain_kind merge_expr_domain
(
    const expr_domain_kind lhsDomain,
    const expr_domain_kind rhsDomain
)
{
    return
        lhsDomain == expr_domain_kind::agnostic ? rhsDomain
      : rhsDomain == expr_domain_kind::agnostic ? lhsDomain
      : lhsDomain == rhsDomain ? lhsDomain
      : expr_domain_kind::incompatible;
}


/**
 * \brief Compile-time domain carried by an expression node.
 *
 * This is the enum counterpart of \c expr_result_t: it extracts one key piece
 * of metadata from an expression type so other templates can reason about it.
 */
template<class Expr>
inline constexpr expr_domain_kind expr_domain_v = std::decay_t<Expr>::domain;


/**
 * \brief Check whether two expression domains may appear in the same tree.
 *
 * For non-expression types the default answer is true, because literals will
 * later be promoted to \c uniform_expr and treated as domain-agnostic.
 *
 * Once both sides are actual expression nodes, compatibility is delegated to
 * \c merge_expr_domain().
 */
template<class LExpr, class RExpr, bool = is_expression_v<LExpr> && is_expression_v<RExpr>>
struct is_expr_domain_compatible
:
    std::true_type
{};


template<class LExpr, class RExpr>
struct is_expr_domain_compatible<LExpr, RExpr, true>
:
    std::integral_constant
    <
        bool,
        merge_expr_domain(expr_domain_v<LExpr>, expr_domain_v<RExpr>)
            != expr_domain_kind::incompatible
    >
{};


template<class LExpr, class RExpr>
inline constexpr bool is_expr_domain_compatible_v =
    is_expr_domain_compatible<LExpr, RExpr>::value;

} // End namespace dg
} // End namespace Foam

#endif
