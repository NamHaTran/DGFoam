/**
 * \file
 * \brief Core traits and domain metadata for the \c Foam::dg expression system.
 */

#ifndef Foam_dgExprTraits_H
#define Foam_dgExprTraits_H

namespace Foam
{
namespace dg
{

/**
 * \brief Detect whether a type participates in the DG expression system.
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
 */
template<class T>
using expr_result_t = typename std::decay_t<T>::value_type;


/**
 * \brief Detect values that may be promoted to a uniform expression leaf.
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
 */
template<class T>
struct is_expr_operand
:
    std::integral_constant<bool, is_expression_v<T> || is_uniform_value_v<T>>
{};


template<class T>
inline constexpr bool is_expr_operand_v = is_expr_operand<T>::value;


/// Sentinel used when an expression leaf cannot report its extent directly.
static constexpr label unknown_expr_extent = -1;


/**
 * \brief Logical domain covered by an expression node.
 *
 * This keeps cell-only, face-only and full Gauss-field expressions from being
 * mixed accidentally.
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


/// Convenience accessor for whether an expression can evaluate cell values.
template<class Expr>
inline constexpr bool expr_supports_cell_v = std::decay_t<Expr>::supports_cell;


/// Convenience accessor for whether an expression can evaluate face values.
template<class Expr>
inline constexpr bool expr_supports_face_v = std::decay_t<Expr>::supports_face;


/**
 * \brief Merge two expression domains into their compatible common domain.
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


/// Compile-time domain carried by an expression node.
template<class Expr>
inline constexpr expr_domain_kind expr_domain_v = std::decay_t<Expr>::domain;


/**
 * \brief Check whether two expression domains may appear in the same tree.
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
