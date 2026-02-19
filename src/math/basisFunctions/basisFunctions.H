/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::math::basisFunctions

Description
    Utility functions for evaluating DG basis functions and their
    derivatives on reference cells and faces.

SourceFiles
    basisFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef basisFunctions_H
#define basisFunctions_H

#include <boost/math/special_functions/jacobi.hpp>
#include <boost/math/tools/config.hpp>
#include "dgCellType.H"
#include "dgFacePosition.H"
#include "scalar.H"
#include "label.H"
#include "List.H"
#include "error.H"
#include "refCoordTransforms.H"

namespace Foam
{

/**
 * \brief Return number of basis functions
 *
 * \param pOrder Polynomial order
 * \param type   DG cell type
 * \return       Number of basis functions
 */
inline label getNumBasis(const label pOrder, const dgCellType type)
{
    label count = 0;

    switch (type)
    {
        // --------------------------------------------------
        // Hexahedron
        // Constraint: i + j + k <= r
        // --------------------------------------------------
        case dgCellType::HEX:
        {
            for (label i = 0; i <= pOrder; ++i)
                for (label j = 0; j <= pOrder - i; ++j)
                    for (label k = 0; k <= pOrder - i - j; ++k)
                        ++count;

            return count;
        }

        // --------------------------------------------------
        // Prism (wedge)
        // Constraint: i <= r,  j + k <= r
        // --------------------------------------------------
        case dgCellType::PRISM:
        {
            for (label i = 0; i <= pOrder; ++i)
                for (label j = 0; j <= pOrder; ++j)
                    for (label k = 0; k <= pOrder - j; ++k)
                        ++count;

            return count;
        }

        // --------------------------------------------------
        // Pyramid
        // Constraint: 0 <= i,j <= r - k,  k <= r
        // --------------------------------------------------
        case dgCellType::PYRAMID:
        {
            for (label i = 0; i <= pOrder; ++i)
                for (label j = 0; j <= pOrder; ++j)
                    for (label k = 0; k <= pOrder - j; ++k)
                        ++count;

            return count;
        }

        // --------------------------------------------------
        // Tetrahedron
        // Constraint: i + j + k <= r
        // --------------------------------------------------
        case dgCellType::TET:
        {
            for (label i = 0; i <= pOrder; ++i)
                for (label j = 0; j <= pOrder - i; ++j)
                    for (label k = 0; k <= pOrder - i - j; ++k)
                        ++count;

            return count;
        }

        case dgCellType::INVALID:
        default:
        {
            FatalErrorInFunction
                << "Unsupported or invalid cell type for basis function: "
                << type << nl
                << exit(FatalError);

            return -1; // suppress compiler warning
        }
    }
}


/**
 * \brief Evaluate orthogonal interior basis functions for hexahedral cells
 *
 * This routine evaluates an orthogonal modal basis defined on the
 * reference hexahedron. The basis functions are constructed as
 * tensor products of Legendre polynomials.
 *
 * The modal basis functions are defined as:
 *
 * \f[
 *   \psi_{i,j,k}(x,y,z)
 *   = P_i(x)\, P_j(y)\, P_k(z),
 *   \qquad i + j + k \le p
 * \f]
 *
 * where \f$P_n\f$ denotes the Legendre polynomial of degree \f$n\f$.
 *
 * The reference domain of the basis is the unit hexahedron:
 *
 * \f[
 *   (x,y,z) \in [-1,1]^3
 * \f]
 *
 * The basis is orthogonal with respect to the \f$L^2\f$ inner product
 * over the reference domain.
 *
 * The routine evaluates:
 *  - the values of all modal basis functions at a given point,
 *  - and their first derivatives with respect to each coordinate.
 *
 * \param etaPt        Reference coordinates in the basis domain
 * \param pOrder       Polynomial order of the approximation
 * \param basis        Values of the basis functions at the evaluation point
 * \param dBasis_dx    Derivatives with respect to x
 * \param dBasis_dy    Derivatives with respect to y
 * \param dBasis_dz    Derivatives with respect to z
 */
inline void computeHexBasisAndDerivatives
(
    const vector etaPt,
    const label  pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_dx,
    List<scalar>& dBasis_dy,
    List<scalar>& dBasis_dz
)
{
    // Reference coordinates transformation:
    // vector X = xiToXHex(etaToXiHex(etaPt));
    scalar x = etaPt.x();
    scalar y = etaPt.y();
    scalar z = etaPt.z();


    // --------------------------------------------------
    // Orthogonal Legendre basis on reference hexahedron
    // Basis definition:
    //   psi_{i,j,k} = P_i(x) P_j(y) P_k(z)
    //   with i + j + k <= pOrder
    // --------------------------------------------------

    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    // Number of basis functions for total-degree HEX
    const label nBasis = getNumBasis(pOrder, dgCellType::HEX);

    basis.setSize(nBasis);
    dBasis_dx.setSize(nBasis);
    dBasis_dy.setSize(nBasis);
    dBasis_dz.setSize(nBasis);

    label idx = 0;

    // Loop over total degree
    for (label i = 0; i <= pOrder; ++i)
    {
        const scalar Pi  = jacobi(static_cast<unsigned>(i), 0.0, 0.0, x);
        const scalar dPi = jacobi_prime(static_cast<unsigned>(i), 0.0, 0.0, x);

        for (label j = 0; j <= pOrder - i; ++j)
        {
            const scalar Pj  = jacobi(static_cast<unsigned>(j), 0.0, 0.0, y);
            const scalar dPj = jacobi_prime(static_cast<unsigned>(j), 0.0, 0.0, y);

            for (label k = 0; k <= pOrder - i - j; ++k)
            {
                const scalar Pk  = jacobi(static_cast<unsigned>(k), 0.0, 0.0, z);
                const scalar dPk = jacobi_prime(static_cast<unsigned>(k), 0.0, 0.0, z);

                // Basis value
                basis[idx] = Pi * Pj * Pk;

                // Derivatives w.r.t reference coordinates
                dBasis_dx[idx] = dPi * Pj  * Pk;
                dBasis_dy[idx] = Pi  * dPj * Pk;
                dBasis_dz[idx] = Pi  * Pj  * dPk;

                ++idx;
            }
        }
    }
}


/**
 * \brief Evaluate orthogonal modal basis functions for prism (wedge) cells
 *
 * This routine evaluates a set of hierarchical orthogonal polynomial
 * basis functions defined on the reference prism (wedge) element.
 * The basis is constructed such that it is orthogonal with respect
 * to the L² inner product on the reference prism when combined with
 * the appropriate Duffy-type quadrature.
 *
 * The orthogonal basis functions are defined as:
 *
 * \f[
 *   \phi_{i,j,k}(x,y,z)
 *   =
 *   P_i^{0,0}(x)\;
 *   P_j^{0,0}(y)\;
 *   \left(\frac{1 - z}{2}\right)^j
 *   P_k^{2j+1,0}(z),
 * \f]
 *
 * with index constraints
 *
 * \f[
 *   i \le p, \qquad
 *   j \le p, \qquad
 *   j + k \le p.
 * \f]
 *
 * Here, \f$P_n^{\alpha,\beta}\f$ denotes the Jacobi polynomial of degree
 * \f$n\f$ with parameters \f$(\alpha,\beta)\f$.
 *
 * The reference domain is the tensor-product prism:
 *
 * \f[
 *   (x,y,z) \in [-1,1] \times [-1,1] \times [-1,1],
 * \f]
 *
 * where the triangular structure in the (y,z) directions is enforced
 * through the factor \f$\left(\frac{1 - z}{2}\right)^j\f$ and the
 * constraint \f$j + k \le p\f$.
 *
 * The resulting basis forms a hierarchical polynomial space suitable
 * for discontinuous Galerkin discretizations on prism elements.
 *
 * \param etaPt        Reference coordinates (x,y,z) in [-1,1]^3
 * \param pOrder       Polynomial order of the approximation
 * \param basis        Values of the basis functions at etaPt
 * \param dBasis_dx    Derivatives with respect to x
 * \param dBasis_dy    Derivatives with respect to y
 * \param dBasis_dz    Derivatives with respect to z
 */
inline void computePrismBasisAndDerivatives
(
    const vector etaPt,
    const label  pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_dx,
    List<scalar>& dBasis_dy,
    List<scalar>& dBasis_dz
)
{
    const scalar x = etaPt.x();
    const scalar y = etaPt.y();
    const scalar z = etaPt.z();

    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    const label nBasis = getNumBasis(pOrder, dgCellType::PRISM);

    basis.setSize(nBasis);
    dBasis_dx.setSize(nBasis);
    dBasis_dy.setSize(nBasis);
    dBasis_dz.setSize(nBasis);

    label idx = 0;

    const scalar half = 0.5;
    const scalar oneMinusZ = 1.0 - z;

    for (label i = 0; i <= pOrder; ++i)
    {
        const scalar Pi  = jacobi(static_cast<unsigned>(i), 0.0, 0.0, x);
        const scalar dPi = jacobi_prime(static_cast<unsigned>(i), 0.0, 0.0, x);

        for (label j = 0; j <= pOrder; ++j)
        {
            const scalar Pj  = jacobi(static_cast<unsigned>(j), 0.0, 0.0, y);
            const scalar dPj = jacobi_prime(static_cast<unsigned>(j), 0.0, 0.0, y);

            const scalar S =
                (j == 0 ? 1.0 : std::pow(half*oneMinusZ, j));

            const scalar dS_dz =
                (j == 0 ? 0.0 :
                 -j * std::pow(half*oneMinusZ, j-1) * half);

            for (label k = 0; k <= pOrder - j; ++k)
            {
                const scalar alpha = 2.0*j + 1.0;

                const scalar Pk  =
                    jacobi(static_cast<unsigned>(k), alpha, 0.0, z);

                const scalar dPk =
                    jacobi_prime(static_cast<unsigned>(k), alpha, 0.0, z);

                // =============================
                // Basis
                // =============================
                basis[idx] = Pi * Pj * S * Pk;

                // =============================
                // d/dx
                // =============================
                dBasis_dx[idx] = dPi * Pj * S * Pk;

                // =============================
                // d/dy
                // =============================
                dBasis_dy[idx] = Pi * dPj * S * Pk;

                // =============================
                // d/dz
                // =============================
                dBasis_dz[idx] =
                    Pi * Pj * (S * dPk + Pk * dS_dz);

                ++idx;
            }
        }
    }
}


/**
 * \brief Evaluate orthogonal interior basis functions for pyramid cells
 *
 * This routine evaluates an orthogonal modal basis defined on the
 * reference pyramid. The polynomial space is identical to that of the
 * classical high-order pyramid element, but the basis functions are
 * constructed to be orthogonal with respect to the L^2 inner product.
 *
 * The orthogonal modal basis functions are defined as:
 *
 * \f[
 *   \psi_{i,j,k}(\eta_1,\eta_2,\eta_3)
 *   =
 *   P_i^{0,0}(\eta_1)\,
 *   P_j^{0,0}(\eta_2)\,
 *   \left(\frac{1-\eta_3}{2}\right)^{\max(i,j)}\,
 *   P_k^{2\max(i,j)+2,\,0}(\eta_3),
 * \f]
 *
 * with index ranges:
 *
 * \f[
 *   0 \le i \le p,\quad
 *   0 \le j \le p,\quad
 *   0 \le k \le p - j,
 * \f]
 *
 * where \f$P_n^{\alpha,\beta}\f$ denotes the Jacobi polynomial of degree
 * \f$n\f$ with parameters \f$(\alpha,\beta)\f$.
 *
 * The reference domain is the cube:
 *
 * \f[
 *   (\eta_1,\eta_2,\eta_3) \in [-1,1]^3,
 * \f]
 *
 * together with the pyramid collapse structure implied by the
 * factor \f$((1-\eta_3)/2)^{\max(i,j)}\f$, which ensures regular
 * behavior at the apex.
 *
 * The routine evaluates:
 *  - the values of all modal basis functions at a given point,
 *  - and their first derivatives with respect to each reference coordinate.
 *
 * \param etaPt        Reference coordinates (\f$\eta_1,\eta_2,\eta_3\f$)
 * \param pOrder       Polynomial order of the approximation
 * \param basis        Values of the basis functions at the evaluation point
 * \param dBasis_deta1 Derivatives with respect to \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives with respect to \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives with respect to \f$\eta_3\f$
 */
inline void computePyramidBasisAndDerivatives
(
    const vector etaPt,     // (eta1, eta2, eta3) in [-1,1]^3
    const label  pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    const scalar eta1 = etaPt.x();
    const scalar eta2 = etaPt.y();
    const scalar eta3 = etaPt.z();

    const label nBasis = getNumBasis(pOrder, dgCellType::PYRAMID);

    basis.setSize(nBasis);
    dBasis_deta1.setSize(nBasis);
    dBasis_deta2.setSize(nBasis);
    dBasis_deta3.setSize(nBasis);

    label idx = 0;

    const scalar oneMinusEta3 = 1.0 - eta3;
    const scalar halfFactor   = 0.5;

    for (label i = 0; i <= pOrder; ++i)
    {
        const scalar Pi  = jacobi(static_cast<unsigned>(i), 0.0, 0.0, eta1);
        const scalar dPi = jacobi_prime(static_cast<unsigned>(i), 0.0, 0.0, eta1);

        for (label j = 0; j <= pOrder; ++j)
        {
            const scalar Pj  = jacobi(static_cast<unsigned>(j), 0.0, 0.0, eta2);
            const scalar dPj = jacobi_prime(static_cast<unsigned>(j), 0.0, 0.0, eta2);

            const label m = max(i, j);

            // ((1 - eta3)/2)^m
            const scalar baseTerm =
                (m == 0)
                ? 1.0
                : std::pow(oneMinusEta3 * halfFactor, m);

            const scalar dBaseTerm_deta3 =
                (m == 0)
                ? 0.0
                : -m * halfFactor *
                  std::pow(oneMinusEta3 * halfFactor, m - 1);

            // --------------------------------------------------
            // Corrected constraint:
            // k ≤ pOrder - j
            // --------------------------------------------------
            for (label k = 0; k <= pOrder - j; ++k)
            {
                const scalar alpha = 2.0 * m + 2.0;

                const scalar Pk =
                    jacobi(static_cast<unsigned>(k), alpha, 0.0, eta3);

                const scalar dPk =
                    jacobi_prime(static_cast<unsigned>(k), alpha, 0.0, eta3);

                // --------------------------------------------------
                // Basis value
                // --------------------------------------------------
                basis[idx] = Pi * Pj * baseTerm * Pk;

                // --------------------------------------------------
                // ∂/∂eta1
                // --------------------------------------------------
                dBasis_deta1[idx] =
                    dPi * Pj * baseTerm * Pk;

                // --------------------------------------------------
                // ∂/∂eta2
                // --------------------------------------------------
                dBasis_deta2[idx] =
                    Pi * dPj * baseTerm * Pk;

                // --------------------------------------------------
                // ∂/∂eta3
                // --------------------------------------------------
                dBasis_deta3[idx] =
                      Pi * Pj * dBaseTerm_deta3 * Pk
                    + Pi * Pj * baseTerm * dPk;

                ++idx;
            }
        }
    }
}

/**
 * \brief Evaluate orthogonal interior basis functions for tetrahedral cells
 *
 * This routine evaluates an orthogonal modal basis defined on the
 * reference tetrahedron. The polynomial space corresponds to the
 * classical tetrahedral finite element space of total degree less
 * than or equal to p, but the basis is constructed to be orthogonal
 * with respect to the L^2 inner product.
 *
 * The orthogonal modal basis functions are defined as:
 *
 * \f[
 *   \psi_{i,j,k}(\eta_1,\eta_2,\eta_3)
 *   =
 *   P_i^{0,0}(\eta_1)
 *   \left(\frac{1-\eta_2}{2}\right)^i
 *   P_j^{2i+1,0}(\eta_2)
 *   \left(\frac{1-\eta_3}{2}\right)^j
 *   P_k^{2i+2j+2,0}(\eta_3),
 *   \qquad i + j + k \le p,
 * \f]
 *
 * where \f$P_n^{\alpha,\beta}\f$ denotes the Jacobi polynomial of degree
 * \f$n\f$ with parameters \f$(\alpha,\beta)\f$.
 *
 * The reference domain is the cube:
 *
 * \f[
 *   (\eta_1,\eta_2,\eta_3) \in [-1,1]^3,
 * \f]
 *
 * together with the tetrahedral collapse structure induced by the
 * hierarchical factors \f$((1-\eta_2)/2)^i\f$ and
 * \f$((1-\eta_3)/2)^j\f$, which ensure polynomial completeness and
 * regular behavior toward the tetrahedral apex.
 *
 * The routine evaluates:
 *  - the values of all modal basis functions at a given point,
 *  - and their first derivatives with respect to each reference coordinate.
 *
 * \param etaPt          Reference coordinates (\f$\eta_1,\eta_2,\eta_3\f$)
 * \param pOrder         Polynomial order of the approximation
 * \param basis          Values of the basis functions at the evaluation point
 * \param dBasis_deta1   Derivatives with respect to \f$\eta_1\f$
 * \param dBasis_deta2   Derivatives with respect to \f$\eta_2\f$
 * \param dBasis_deta3   Derivatives with respect to \f$\eta_3\f$
 */
inline void computeTetBasisAndDerivatives
(
    const vector etaPt,   // (eta1, eta2, eta3) in [-1,1]^3
    const label  pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    const scalar eta1 = etaPt.x();
    const scalar eta2 = etaPt.y();
    const scalar eta3 = etaPt.z();

    const label nBasis = getNumBasis(pOrder, dgCellType::TET);

    basis.setSize(nBasis);
    dBasis_deta1.setSize(nBasis);
    dBasis_deta2.setSize(nBasis);
    dBasis_deta3.setSize(nBasis);

    label idx = 0;

    const scalar half = 0.5;

    for (label i = 0; i <= pOrder; ++i)
    {
        const scalar Pi  = jacobi(static_cast<unsigned>(i), 0.0, 0.0, eta1);
        const scalar dPi = jacobi_prime(static_cast<unsigned>(i), 0.0, 0.0, eta1);

        const scalar oneMinusEta2 = 1.0 - eta2;

        // ((1 - eta2)/2)^i
        const scalar powEta2 =
            (i == 0) ? 1.0 :
            std::pow(oneMinusEta2 * half, i);

        const scalar dPowEta2_deta2 =
            (i == 0) ? 0.0 :
            -i * half *
            std::pow(oneMinusEta2 * half, i - 1);

        for (label j = 0; j <= pOrder - i; ++j)
        {
            const scalar alpha_j = 2.0 * i + 1.0;

            const scalar Pj =
                jacobi(static_cast<unsigned>(j), alpha_j, 0.0, eta2);

            const scalar dPj =
                jacobi_prime(static_cast<unsigned>(j), alpha_j, 0.0, eta2);

            const scalar oneMinusEta3 = 1.0 - eta3;

            // ((1 - eta3)/2)^j
            const scalar powEta3 =
                (j == 0) ? 1.0 :
                std::pow(oneMinusEta3 * half, j);

            const scalar dPowEta3_deta3 =
                (j == 0) ? 0.0 :
                -j * half *
                std::pow(oneMinusEta3 * half, j - 1);

            for (label k = 0; k <= pOrder - i - j; ++k)
            {
                const scalar alpha_k =
                    2.0 * i + 2.0 * j + 2.0;

                const scalar Pk =
                    jacobi(static_cast<unsigned>(k), alpha_k, 0.0, eta3);

                const scalar dPk =
                    jacobi_prime(static_cast<unsigned>(k), alpha_k, 0.0, eta3);

                // --------------------------------------------------
                // Basis value
                // --------------------------------------------------
                basis[idx] =
                    Pi * powEta2 *
                    Pj * powEta3 *
                    Pk;

                // --------------------------------------------------
                // ∂/∂eta1
                // --------------------------------------------------
                dBasis_deta1[idx] =
                    dPi * powEta2 *
                    Pj * powEta3 *
                    Pk;

                // --------------------------------------------------
                // ∂/∂eta2
                // --------------------------------------------------
                dBasis_deta2[idx] =
                      Pi * dPowEta2_deta2 *
                        Pj * powEta3 * Pk
                    + Pi * powEta2 *
                        dPj * powEta3 * Pk;

                // --------------------------------------------------
                // ∂/∂eta3
                // --------------------------------------------------
                dBasis_deta3[idx] =
                      Pi * powEta2 *
                        Pj * dPowEta3_deta3 * Pk
                    + Pi * powEta2 *
                        Pj * powEta3 * dPk;

                ++idx;
            }
        }
    }
}

} // namespace Foam

#endif

// ************************************************************************* //
