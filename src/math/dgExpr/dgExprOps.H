/**
 * \file
 * \brief Stateless operator functors used by DG expression nodes.
 */

#ifndef Foam_dgExprOps_H
#define Foam_dgExprOps_H

namespace Foam
{
namespace dg
{

/// Binary addition.
struct add_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs + rhs;
    }
};


/// Binary subtraction.
struct sub_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs - rhs;
    }
};


/// Binary multiplication.
struct mul_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs * rhs;
    }
};


/// Binary division.
struct div_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs / rhs;
    }
};


/// Vector/tensor cross-like product using \c operator^.
struct cross_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs ^ rhs;
    }
};


/// Dot product using \c operator&.
struct dot_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs & rhs;
    }
};


/// Double-dot contraction using \c operator&&.
struct double_dot_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return lhs && rhs;
    }
};


/// Power operation dispatched to \c Foam::pow.
struct pow_op
{
    template<class L, class R>
    static inline auto apply(const L& lhs, const R& rhs)
    {
        return ::Foam::pow(lhs, rhs);
    }
};


/// Unary negation.
struct neg_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return -value;
    }
};


/// Magnitude.
struct mag_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::mag(value);
    }
};


/// Squared magnitude.
struct mag_sqr_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::magSqr(value);
    }
};


/// Square root.
struct sqrt_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::sqrt(value);
    }
};


/// Deviatoric part.
struct dev_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::dev(value);
    }
};


/// Symmetric part.
struct symm_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::symm(value);
    }
};


/// Skew-symmetric part.
struct skew_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::skew(value);
    }
};


/// Transpose.
struct transpose_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::T(value);
    }
};


/// Trace.
struct trace_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::tr(value);
    }
};


/// Determinant.
struct det_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::det(value);
    }
};


/// Inverse.
struct inv_op
{
    template<class T>
    static inline auto apply(const T& value)
    {
        return ::Foam::inv(value);
    }
};

} // End namespace dg
} // End namespace Foam

#endif
