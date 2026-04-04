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
 *
 * From a learning perspective, leaves are the easiest place to start because
 * they reveal the access contract expected by the rest of the tree:
 *
 * - report whether cell and/or face access is supported;
 * - report a size hint when one is known;
 * - provide pointwise getters such as \c cellValue() and
 *   \c faceMinusValue().
 *
 * Internal nodes do not care whether a child is a field reference, a literal
 * or another composed expression. They only care that the child answers this
 * common API.
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
 *
 * This class is the canonical "enter the lazy world" adapter for a full
 * \c GaussField. If user code writes
 *
 *     auto e = dg::expr(U);
 *
 * then \c e is typically a \c gauss_field_ref<Type>.
 *
 * Important mental model:
 *
 * - the actual numerical data still lives in \c field_;
 * - this object only stores how to read that data;
 * - later recipe nodes call \c cellValue() / \c faceValue() on this leaf
 *   whenever they need one quadrature-point sample.
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

    /**
     * \brief Wrap an existing \c GaussField by const reference.
     *
     * No copy of the DG storage is made here.
     */
    explicit gauss_field_ref(const GaussField<Type>& field)
    :
        field_(field)
    {}

    /// \brief Number of cell-interior Gauss points available from \c field_.
    inline label cellSizeHint() const
    {
        return field_.cellField().size();
    }

    /// \brief Total number of face Gauss points available from \c field_.
    inline label faceSizeHint() const
    {
        return field_.faceField().nGauss();
    }

    /// \brief Read one interior Gauss-point value from the wrapped field.
    inline const Type& cellValue(const label gpI) const
    {
        return field_.cellField()[gpI];
    }

    /// \brief Read the interior-side trace at face Gauss-point \c gpI.
    inline const Type& faceMinusValue(const label gpI) const
    {
        return field_.faceField().minusValue(gpI);
    }

    /// \brief Read the exterior-side trace at face Gauss-point \c gpI.
    inline const Type& facePlusValue(const label gpI) const
    {
        return field_.faceField().plusValue(gpI);
    }

    /**
     * \brief Read the interior-side trace on one specific local face.
     *
     * This accessor is useful for APIs that work face-by-face instead of over
     * the flat face-Gauss storage.
     */
    inline const Type& faceMinusValueOnFace
    (
        const label faceI,
        const label gpI
    ) const
    {
        return field_.faceField().minusValueOnFace(faceI, gpI);
    }

    /**
     * \brief Read the exterior-side trace on one specific local face.
     */
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
 *
 * This is the right leaf when user code enters the expression system through
 * \c someGaussField.cellField().
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
 *
 * This is typically produced by \c dg::expr(someGaussField.faceField()).
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
 *
 * Its size hints are intentionally unknown. A literal gets its runtime extent
 * only after it is paired with a field-backed expression.
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
 *
 * This is the ownership-preserving counterpart of \c gauss_field_ref.
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
 * This alias is one of the key bridges between the user-facing API and the
 * internal recipe representation.
 *
 * Rules:
 *
 * - existing expressions stay unchanged;
 * - plain literals become \c uniform_expr.
 *
 * Example:
 *
 *     expression_node_t<gauss_field_ref<scalar>> -> gauss_field_ref<scalar>
 *     expression_node_t<double> -> uniform_expr<double>
 *
 * The operator builders use this alias so they can accept mixed inputs such as
 * \c dg::expr(p) + 1.0 while still storing a uniform internal node type.
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
 *
 * This overload preserves the lazy node that already exists.
 */
template<class T>
inline auto make_expression_node(const T& value)
    -> std::enable_if_t<is_expression_v<T>, std::decay_t<T>>
{
    return std::decay_t<T>(value);
}


/**
 * \brief Promote a literal value to a uniform expression leaf.
 *
 * This overload is what lets user algebra mix constants with expressions
 * without having to spell \c dg::uniform(...) explicitly.
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
