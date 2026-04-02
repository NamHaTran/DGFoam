/**
 * \file
 * \brief Leaf-node adapters that expose DG field storage to the expression
 *        template layer without copying.
 *
 * Leaves are the bridge between "real storage" and "lazy expression nodes".
 *
 * A practical way to read this file:
 *
 * - leaf classes answer "where does one value come from?"
 * - they do not decide how values are combined
 * - they expose a uniform access API used by the rest of the expression tree
 *
 * Most leaves keep only a reference to an existing field, so building an
 * expression does not duplicate DG data. The exception is literals and
 * temporaries, which must store a value or a \c tmp handle to stay valid until
 * the expression is eventually evaluated.
 */

#ifndef Foam_dgExprLeaves_H
#define Foam_dgExprLeaves_H

namespace Foam
{
namespace dg
{

/**
 * \brief Non-owning expression leaf for a full \c GaussField.
 *
 * Typical construction:
 *
 *     auto e = dg::expr(U);
 *
 * This object does not own any quadrature values. It only stores a
 * \c const\ GaussField<Type>& and forwards all reads to the referenced field.
 *
 * That means:
 *
 * - creating the expression is cheap;
 * - no cell/face arrays are allocated here;
 * - the original field must outlive the expression tree.
 */
template<class Type>
class gauss_field_ref
{
    const GaussField<Type>& field_;

public:

    using expression_tag = void;
    using value_type = Type;
    static constexpr bool supports_cell = true;
    static constexpr bool supports_face = true;
    static constexpr expr_domain_kind domain = expr_domain_kind::gauss;

    explicit gauss_field_ref(const GaussField<Type>& field)
    :
        field_(field)
    {}

    inline label cellSizeHint() const
    {
        return field_.cellField().size();
    }

    inline label faceSizeHint() const
    {
        return field_.faceField().nGauss();
    }

    inline const Type& cellValue(const label gpI) const
    {
        return field_.cellField()[gpI];
    }

    inline const Type& faceMinusValue(const label gpI) const
    {
        return field_.faceField().minusValue(gpI);
    }

    inline const Type& facePlusValue(const label gpI) const
    {
        return field_.faceField().plusValue(gpI);
    }

    inline const Type& faceMinusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.faceField().minusValueOnFace(faceI, gpI);
    }

    inline const Type& facePlusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.faceField().plusValueOnFace(faceI, gpI);
    }
};


/**
 * \brief Non-owning expression leaf for cell-interior Gauss values only.
 *
 * The compile-time flags advertise that only cell values are available.
 * Operations that require face access can therefore be rejected or excluded
 * without runtime checks.
 */
template<class Type>
class cell_gauss_field_ref
{
    const cellGaussField<Type>& field_;

public:

    using expression_tag = void;
    using value_type = Type;
    static constexpr bool supports_cell = true;
    static constexpr bool supports_face = false;
    static constexpr expr_domain_kind domain = expr_domain_kind::cell;

    explicit cell_gauss_field_ref(const cellGaussField<Type>& field)
    :
        field_(field)
    {}

    inline label cellSizeHint() const
    {
        return field_.size();
    }

    inline label faceSizeHint() const
    {
        return unknown_expr_extent;
    }

    inline const Type& cellValue(const label gpI) const
    {
        return field_[gpI];
    }
};


/**
 * \brief Non-owning expression leaf for face Gauss values only.
 *
 * This is the face-side equivalent of \c cell_gauss_field_ref. It lets the
 * tree work with face data while clearly stating that no cell values exist.
 */
template<class Type>
class face_gauss_field_ref
{
    const faceGaussField<Type>& field_;

public:

    using expression_tag = void;
    using value_type = Type;
    static constexpr bool supports_cell = false;
    static constexpr bool supports_face = true;
    static constexpr expr_domain_kind domain = expr_domain_kind::face;

    explicit face_gauss_field_ref(const faceGaussField<Type>& field)
    :
        field_(field)
    {}

    inline label cellSizeHint() const
    {
        return unknown_expr_extent;
    }

    inline label faceSizeHint() const
    {
        return field_.nGauss();
    }

    inline const Type& faceMinusValue(const label gpI) const
    {
        return field_.minusValue(gpI);
    }

    inline const Type& facePlusValue(const label gpI) const
    {
        return field_.plusValue(gpI);
    }

    inline const Type& faceMinusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.minusValueOnFace(faceI, gpI);
    }

    inline const Type& facePlusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.plusValueOnFace(faceI, gpI);
    }
};


/**
 * \brief Expression leaf that repeats a single literal value everywhere.
 *
 * Example:
 *
 *     auto e = dg::expr(p) + 1.0;
 *
 * The \c 1.0 is normalized into a \c uniform_expr<scalar>. Unlike field-based
 * leaves, this class stores its own copy because a literal has no external
 * storage to reference.
 */
template<class Type>
class uniform_expr
{
    Type value_;

public:

    using expression_tag = void;
    using value_type = Type;
    static constexpr bool supports_cell = true;
    static constexpr bool supports_face = true;
    static constexpr expr_domain_kind domain = expr_domain_kind::agnostic;

    explicit uniform_expr(const Type& value)
    :
        value_(value)
    {}

    inline label cellSizeHint() const
    {
        return unknown_expr_extent;
    }

    inline label faceSizeHint() const
    {
        return unknown_expr_extent;
    }

    inline const Type& cellValue(const label) const
    {
        return value_;
    }

    inline const Type& faceMinusValue(const label) const
    {
        return value_;
    }

    inline const Type& facePlusValue(const label) const
    {
        return value_;
    }

    inline const Type& faceMinusValueOnFace(const label, const label) const
    {
        return value_;
    }

    inline const Type& facePlusValueOnFace(const label, const label) const
    {
        return value_;
    }
};


/**
 * \brief Owning expression leaf for a temporary full \c GaussField.
 *
 * The wrapped \c tmp keeps the evaluated storage alive for the lifetime of the
 * expression tree.
 *
 * This is important when an expression is built from a temporary:
 *
 *     auto e = dg::expr(makeTemporaryField());
 *
 * A reference-only leaf would dangle as soon as the temporary was destroyed.
 * By storing the \c tmp object, the leaf keeps the temporary storage alive
 * until the whole lazy expression is consumed.
 */
template<class Type>
class tmp_gauss_field_ref
{
    tmp<GaussField<Type>> field_;

public:

    using expression_tag = void;
    using value_type = Type;
    static constexpr bool supports_cell = true;
    static constexpr bool supports_face = true;
    static constexpr expr_domain_kind domain = expr_domain_kind::gauss;

    explicit tmp_gauss_field_ref(const tmp<GaussField<Type>>& field)
    :
        field_(field)
    {}

    inline label cellSizeHint() const
    {
        return field_.cref().cellField().size();
    }

    inline label faceSizeHint() const
    {
        return field_.cref().faceField().nGauss();
    }

    inline const Type& cellValue(const label gpI) const
    {
        return field_.cref().cellField()[gpI];
    }

    inline const Type& faceMinusValue(const label gpI) const
    {
        return field_.cref().faceField().minusValue(gpI);
    }

    inline const Type& facePlusValue(const label gpI) const
    {
        return field_.cref().faceField().plusValue(gpI);
    }

    inline const Type& faceMinusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.cref().faceField().minusValueOnFace(faceI, gpI);
    }

    inline const Type& facePlusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.cref().faceField().plusValueOnFace(faceI, gpI);
    }
};


/**
 * \brief Owning expression leaf for a temporary \c cellGaussField.
 *
 * Same ownership idea as \c tmp_gauss_field_ref, but restricted to cell data.
 */
template<class Type>
class tmp_cell_gauss_field_ref
{
    tmp<cellGaussField<Type>> field_;

public:

    using expression_tag = void;
    using value_type = Type;
    static constexpr bool supports_cell = true;
    static constexpr bool supports_face = false;
    static constexpr expr_domain_kind domain = expr_domain_kind::cell;

    explicit tmp_cell_gauss_field_ref(const tmp<cellGaussField<Type>>& field)
    :
        field_(field)
    {}

    inline label cellSizeHint() const
    {
        return field_.cref().size();
    }

    inline label faceSizeHint() const
    {
        return unknown_expr_extent;
    }

    inline const Type& cellValue(const label gpI) const
    {
        return field_.cref()[gpI];
    }
};


/**
 * \brief Owning expression leaf for a temporary \c faceGaussField.
 *
 * Same ownership idea as \c tmp_gauss_field_ref, but restricted to face data.
 */
template<class Type>
class tmp_face_gauss_field_ref
{
    tmp<faceGaussField<Type>> field_;

public:

    using expression_tag = void;
    using value_type = Type;
    static constexpr bool supports_cell = false;
    static constexpr bool supports_face = true;
    static constexpr expr_domain_kind domain = expr_domain_kind::face;

    explicit tmp_face_gauss_field_ref(const tmp<faceGaussField<Type>>& field)
    :
        field_(field)
    {}

    inline label cellSizeHint() const
    {
        return unknown_expr_extent;
    }

    inline label faceSizeHint() const
    {
        return field_.cref().nGauss();
    }

    inline const Type& faceMinusValue(const label gpI) const
    {
        return field_.cref().minusValue(gpI);
    }

    inline const Type& facePlusValue(const label gpI) const
    {
        return field_.cref().plusValue(gpI);
    }

    inline const Type& faceMinusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.cref().minusValueOnFace(faceI, gpI);
    }

    inline const Type& facePlusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.cref().plusValueOnFace(faceI, gpI);
    }
};


/**
 * \brief Normalize an operand into the expression node type stored by trees.
 *
 * Existing expressions are kept as-is; plain literals are wrapped in
 * \c uniform_expr.
 */
template<class T>
using expression_node_t = std::conditional_t
<
    is_expression_v<T>,
    std::decay_t<T>,
    uniform_expr<std::decay_t<T>>
>;


/**
 * \brief Forward an existing expression node unchanged.
 */
template<class T>
inline auto make_expression_node(const T& value)
    -> std::enable_if_t<is_expression_v<T>, std::decay_t<T>>
{
    return std::decay_t<T>(value);
}


/**
 * \brief Promote a literal value to a uniform expression leaf.
 */
template<class T>
inline auto make_expression_node(const T& value)
    -> std::enable_if_t
    <
        !is_expression_v<T> && is_uniform_value_v<T>,
        uniform_expr<std::decay_t<T>>
    >
{
    return uniform_expr<std::decay_t<T>>(value);
}

} // End namespace dg
} // End namespace Foam

#endif
