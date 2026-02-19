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

namespace Foam
{

/**
 * \brief Return number of interior basis functions
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
        case dgCellType::HEX:
        {
            for (label p = 1; p < pOrder+1; ++p)
                for (label q = 1; q < pOrder+1; ++q)
                    for (label r = 1; r < pOrder+1; ++r)
                        ++count;
            return count;
        }

        case dgCellType::PRISM:
        {
            for (label p = 1; p < pOrder+1; ++p)
                for (label q = 1; q < pOrder+1; ++q)
                    for (label r = 1; r < pOrder+1; ++r)
                        if (q + r < pOrder+1)
                            ++count;
            return count;
        }

        case dgCellType::PYRAMID:
        {
            for (label p = 1; p < pOrder+1; ++p)
                for (label q = 1; q < pOrder+1; ++q)
                    for (label r = 1; r < pOrder+1; ++r)
                        if (p + q + r < pOrder+1)
                            ++count;
            return count;
        }

        case dgCellType::TET:
        {
            for (label p = 1; p < pOrder+1; ++p)
                for (label q = 1; q < pOrder+1; ++q)
                {
                    if (p + q < pOrder+1)
                    {
                        for (label r = 1; r < pOrder+1; ++r)
                        {
                            if (p + q + r < pOrder+1)
                            ++count;
                        }
                    }
                }
            return count;
        }

        case dgCellType::INVALID:
        default:
        {
            FatalErrorInFunction
                << "Unsupported or invalid cell type for basis function: "
                << type << nl
                << exit(FatalError);
            return -1; // Just to satisfy compiler
        }
    }
}

/**
 * \brief Evaluate orthogonal interior basis functions for hexahedral cells
 *
 * The orthogonal basis functions implemented by this routine are defined
 * on a reference hexahedron whose local coordinate system and orientation
 * correspond to the domain illustrated in the figure below.
 * 
 * \image html refHexBasisDomain.png
 *
 * The input coordinates (\f$\eta_1,\eta_2,\eta_3\f$) passed to this function
 * are assumed to lie in the reference coordinate system of that domain.
 *
 * When this routine is used within DGFoam, the Gauss points are typically
 * generated in the solver reference hexahedron
 * \f$(\eta_1,\eta_2,\eta_3) \in [-1,1]^3\f$.
 * In that case, the Gauss-point coordinates must first be mapped to the
 * reference coordinate system of the basis domain using the appropriate
 * coordinate transformation before calling this function.
 *
 * No geometric mapping or coordinate transformation is performed inside
 * this routine. The caller is responsible for providing correctly mapped
 * reference coordinates.
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$ in the basis domain
 * \param eta2         Reference coordinate \f$\eta_2\f$ in the basis domain
 * \param eta3         Reference coordinate \f$\eta_3\f$ in the basis domain
 * \param pOrder       Polynomial order
 * \param basis        Basis function values at the evaluation point
 * \param dBasis_deta1 Derivatives with respect to \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives with respect to \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives with respect to \f$\eta_3\f$
 */
inline void Foam::computeHexBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label  pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    // --------------------------------------------------
    // Orthogonal basis on reference hexahedron [-1,1]^3
    // Tensor-product Legendre polynomials
    // --------------------------------------------------

    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    // Number of basis functions for HEX cell
    const label nBasis = getNumBasis(pOrder, dgCellType::HEX);

    basis.setSize(nBasis);
    dBasis_deta1.setSize(nBasis);
    dBasis_deta2.setSize(nBasis);
    dBasis_deta3.setSize(nBasis);

    label idx = 0;

    // Polynomial degree runs from 0 to pOrder
    for (label p = 0; p <= pOrder; ++p)
    {
        // Legendre polynomial P_p(eta1)
        const scalar Pp  = jacobi(static_cast<unsigned>(p), 0.0, 0.0, eta1);
        const scalar dPp = jacobi_prime(static_cast<unsigned>(p), 0.0, 0.0, eta1);

        for (label q = 0; q <= pOrder; ++q)
        {
            // Legendre polynomial P_q(eta2)
            const scalar Pq  = jacobi(static_cast<unsigned>(q), 0.0, 0.0, eta2);
            const scalar dPq = jacobi_prime(static_cast<unsigned>(q), 0.0, 0.0, eta2);

            for (label r = 0; r <= pOrder; ++r)
            {
                // Legendre polynomial P_r(eta3)
                const scalar Pr  = jacobi(static_cast<unsigned>(r), 0.0, 0.0, eta3);
                const scalar dPr = jacobi_prime(static_cast<unsigned>(r), 0.0, 0.0, eta3);

                // Basis value
                basis[idx] = Pp * Pq * Pr;

                // Derivatives in reference coordinates
                dBasis_deta1[idx] = dPp * Pq  * Pr;
                dBasis_deta2[idx] = Pp  * dPq * Pr;
                dBasis_deta3[idx] = Pp  * Pq  * dPr;

                ++idx;
            }
        }
    }
}

/**
 * \brief Evaluate interior basis functions for prism cells
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$
 * \param eta2         Reference coordinate \f$\eta_2\f$
 * \param eta3         Reference coordinate \f$\eta_3\f$
 * \param pOrder       Polynomial order
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives w.r.t. \f$\eta_3\f$
 */
inline void computeInteriorPrismBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    const label nBasis = getNumBasis(pOrder, dgCellType::PRISM);
    basis.setSize(nBasis);
    dBasis_deta1.setSize(nBasis);
    dBasis_deta2.setSize(nBasis);
    dBasis_deta3.setSize(nBasis);

    label idx = 0;

    label polyDegree = pOrder+1;

    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp     = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp    = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);

        for (label q = 1; q < polyDegree; ++q)
        {
            scalar Pq     = jacobi(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);
            scalar dPq    = jacobi_prime(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);

            const scalar oneMinusEta3 = 1.0 - eta3;
            const scalar oneMinusEta3PowP = std::pow(oneMinusEta3, p);
            const scalar dOneMinusEta3PowP = (p > 0)
                ? -p * std::pow(oneMinusEta3, p - 1)
                : 0.0;

            const double alpha = 2 * q + 1;

            for (label r = 1; r < polyDegree; ++r)
            {
                if (q + r < polyDegree)
                {
                    scalar Pr    = jacobi(static_cast<unsigned>(r-1), alpha, 1.0, eta3);
                    scalar dPr   = jacobi_prime(static_cast<unsigned>(r-1), alpha, 1.0, eta3);

                    // Basis value
                    basis[idx] = Pp * Pq * oneMinusEta3PowP * Pr;

                    // ∂φ/∂eta1
                    dBasis_deta1[idx] = dPp * Pq * oneMinusEta3PowP * Pr;

                    // ∂φ/∂eta2
                    dBasis_deta2[idx] = Pp * dPq * oneMinusEta3PowP * Pr;

                    // ∂φ/∂eta3
                    dBasis_deta3[idx] = Pp * Pq *
                        (dOneMinusEta3PowP * Pr + oneMinusEta3PowP * dPr);

                    ++idx;
                }
            }
        }
    }
}

/**
 * \brief Evaluate interior basis functions for pyramid cells
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$
 * \param eta2         Reference coordinate \f$\eta_2\f$
 * \param eta3         Reference coordinate \f$\eta_3\f$
 * \param pOrder       Polynomial order
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives w.r.t. \f$\eta_3\f$
 */
inline void computeInteriorPyramidBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    const label nBasis = getNumBasis(pOrder, dgCellType::PYRAMID);
    basis.setSize(nBasis);
    dBasis_deta1.setSize(nBasis);
    dBasis_deta2.setSize(nBasis);
    dBasis_deta3.setSize(nBasis);

    label idx = 0;

    label polyDegree = pOrder+1;

    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);

        for (label q = 1; q < polyDegree; ++q)
        {
            scalar Pq  = jacobi(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);
            scalar dPq = jacobi_prime(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);

            double exponent = p + q;
            scalar oneMinusEta3 = 1.0 - eta3;
            scalar powTerm = std::pow(oneMinusEta3, exponent);
            scalar dPowTerm = (exponent > 0)
                ? -exponent * std::pow(oneMinusEta3, exponent - 1)
                : 0.0;

            for (label r = 1; r < polyDegree; ++r)
            {
                if (p + q + r < polyDegree)
                {
                    double alpha = 2 * p + 2 * q + 1;
                    scalar Pr   = jacobi(static_cast<unsigned>(r-1), alpha, 1.0, eta3);
                    scalar dPr  = jacobi_prime(static_cast<unsigned>(r-1), alpha, 1.0, eta3);

                    // Basis value
                    basis[idx] = Pp * Pq * powTerm * Pr;

                    // ∂φ/∂eta1
                    dBasis_deta1[idx] = dPp * Pq * powTerm * Pr;

                    // ∂φ/∂eta2
                    dBasis_deta2[idx] = Pp * dPq * powTerm * Pr;

                    // ∂φ/∂eta3
                    dBasis_deta3[idx] = Pp * Pq * (dPowTerm * Pr + powTerm * dPr);

                    ++idx;
                }
            }
        }
    }
}

/**
 * \brief Evaluate interior basis functions for tetrahedral cells
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$
 * \param eta2         Reference coordinate \f$\eta_2\f$
 * \param eta3         Reference coordinate \f$\eta_3\f$
 * \param pOrder       Polynomial order
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives w.r.t. \f$\eta_3\f$
 */
inline void computeInteriorTetBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    const label nBasis = getNumBasis(pOrder, dgCellType::TET);
    basis.setSize(nBasis);
    dBasis_deta1.setSize(nBasis);
    dBasis_deta2.setSize(nBasis);
    dBasis_deta3.setSize(nBasis);

    label idx = 0;

    label polyDegree = pOrder+1;

    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);

        scalar oneMinusEta2 = 1.0 - eta2;
        scalar oneMinusEta2PowP = std::pow(oneMinusEta2, p);
        scalar dOneMinusEta2PowP = (p > 0) ? -p * std::pow(oneMinusEta2, p - 1) : 0.0;

        for (label q = 1; q < polyDegree; ++q)
        {
            if (p + q < polyDegree)
            {
                double alpha_q = 2 * p + 1;
                scalar Pq     = jacobi(static_cast<unsigned>(q-1), alpha_q, 1.0, eta2);
                scalar dPq    = jacobi_prime(static_cast<unsigned>(q-1), alpha_q, 1.0, eta2);

                scalar oneMinusEta3 = 1.0 - eta3;
                double pqSum = p + q;
                scalar oneMinusEta3PowPQ = std::pow(oneMinusEta3, pqSum);
                scalar dOneMinusEta3PowPQ = (pqSum > 0)
                    ? -pqSum * std::pow(oneMinusEta3, pqSum - 1)
                    : 0.0;

                for (label r = 1; r < polyDegree; ++r)
                {
                    if (p + q + r < polyDegree)
                    {
                        double alpha_r = 2 * p + 2 * q + 1;
                        scalar Pr  = jacobi(static_cast<unsigned>(r-1), alpha_r, 1.0, eta3);
                        scalar dPr = jacobi_prime(static_cast<unsigned>(r-1), alpha_r, 1.0, eta3);

                        // Basis
                        basis[idx] = Pp * oneMinusEta2PowP * Pq * oneMinusEta3PowPQ * Pr;

                        // ∂φ/∂eta1
                        dBasis_deta1[idx] = dPp * oneMinusEta2PowP * Pq * oneMinusEta3PowPQ * Pr;

                        // ∂φ/∂eta2
                        dBasis_deta2[idx] = Pp *
                            (dOneMinusEta2PowP * Pq + oneMinusEta2PowP * dPq) *
                            oneMinusEta3PowPQ * Pr;

                        // ∂φ/∂eta3
                        dBasis_deta3[idx] = Pp * oneMinusEta2PowP * Pq *
                            (dOneMinusEta3PowPQ * Pr + oneMinusEta3PowPQ * dPr);

                        ++idx;
                    }
                }
            }
        }
    }
}

/**
 * \brief Evaluate basis functions on triangular faces
 *
 * \param eta1         Reference face coordinate \f$\eta_1\f$
 * \param eta2         Reference face coordinate \f$\eta_2\f$
 * \param pOrder       Polynomial order
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 */
inline void computeInteriorTriBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2
)
{
/*
Note
    The reference triangle is constructed by collapsing a reference
    quadrilateral along the η₁ direction.

    In this construction, the two vertices of the reference quadrilateral
    are collapsed onto a single vertex along the η₁ axis, forming a
    triangular reference element. As a result, the polynomial basis
    is defined using a tensor-product structure with a collapse factor
    in the η₂ direction.

    When using this function, the reference coordinates must be provided
    consistently with this collapse mapping:
        - η₁ represents the coordinate along the collapsed direction
          (the direction in which the quadrilateral is collapsed to form
          the triangle).
        - η₂ represents the coordinate along the remaining, non-collapsed
          direction.

    Providing (η₁, η₂) in a different order will lead to incorrect
    evaluation of the triangular modal basis functions.
*/

    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    label polyDegree = pOrder+1;

    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        
        const scalar oneMinusEta2 = 1.0 - eta2;
        const scalar oneMinusEta2PowP = std::pow(oneMinusEta2, p);
        const scalar dOneMinusEta2PowP = (p > 0)
            ? -p * std::pow(oneMinusEta2, p - 1)
            : 0.0;
        
        const double alpha = 2 * p + 1;

        for (label q = 1; q < polyDegree; ++q)
        {
            if (p + q < polyDegree)
            {
                scalar Pq  = jacobi(static_cast<unsigned>(q-1), alpha, 1.0, eta2);
                scalar dPq = jacobi_prime(static_cast<unsigned>(q-1), alpha, 1.0, eta2);
                // Basis value
                basis.append(Pp * oneMinusEta2PowP * Pq);
                dBasis_deta1.append(dPp * oneMinusEta2PowP * Pq);
                dBasis_deta2.append(Pp *
                    (dOneMinusEta2PowP * Pq + oneMinusEta2PowP * dPq));
            }
        }
    }
}

/**
 * \brief Evaluate basis functions on quadrilateral faces
 *
 * \param eta1         Reference face coordinate \f$\eta_1\f$
 * \param eta2         Reference face coordinate \f$\eta_2\f$
 * \param pOrder       Polynomial order
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 */
inline void computeInteriorQuadBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const label pOrder,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2
)
{
    using boost::math::jacobi;
    using boost::math::jacobi_prime;

    label polyDegree = pOrder+1;
    
    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp    = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);

        for (label q = 1; q < polyDegree; ++q)
        {
            scalar Pq  = jacobi(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);
            scalar dPq    = jacobi_prime(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);
            // Basis value
            basis.append(Pp * Pq);
            dBasis_deta1.append(dPp * Pq);
            dBasis_deta2.append(Pp * dPq);
        }
    }
}

/**
 * \brief Evaluate face basis functions for hex cells
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$
 * \param eta2         Reference coordinate \f$\eta_2\f$
 * \param eta3         Reference coordinate \f$\eta_3\f$
 * \param pOrder       Polynomial order
 * \param pos          Face position in reference hex
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives w.r.t. \f$\eta_3\f$
 */
inline void computeFaceBasisAndDerivativesOfHex
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    const dgFacePosition pos,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    switch (pos)
    {
        case dgFacePosition::ABCD:
            computeInteriorQuadBasisAndDerivatives(eta1, eta2, pOrder, basis, dBasis_deta1, dBasis_deta2);
            dBasis_deta3.setSize(basis.size());
            forAll(dBasis_deta3, i)
                dBasis_deta3[i] = 0.0;
            break;
        case dgFacePosition::EFGH:
            computeInteriorQuadBasisAndDerivatives(eta1, eta2, pOrder, basis, dBasis_deta1, dBasis_deta2);
            dBasis_deta3.setSize(basis.size());
            forAll(dBasis_deta3, i)
                dBasis_deta3[i] = 0.0;
            break;
        case dgFacePosition::ABEF:
            computeInteriorQuadBasisAndDerivatives(eta1, eta3, pOrder, basis, dBasis_deta1, dBasis_deta3);
            dBasis_deta2.setSize(basis.size());
            forAll(dBasis_deta2, i)
                dBasis_deta2[i] = 0.0;
            break;
        case dgFacePosition::CDGH:
            computeInteriorQuadBasisAndDerivatives(eta1, eta3, pOrder, basis, dBasis_deta1, dBasis_deta3);
            dBasis_deta2.setSize(basis.size());
            forAll(dBasis_deta2, i)
                dBasis_deta2[i] = 0.0;
            break;
        case dgFacePosition::ACEG:
            computeInteriorQuadBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            dBasis_deta1.setSize(basis.size());
            forAll(dBasis_deta1, i)
                dBasis_deta1[i] = 0.0;
            break;
        case dgFacePosition::BDFH:
            computeInteriorQuadBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            dBasis_deta1.setSize(basis.size());
            forAll(dBasis_deta1, i)
                dBasis_deta1[i] = 0.0;
            break;
        default:
            FatalErrorInFunction
                << "Invalid dgFacePosition enum value."
                << abort(FatalError);
    }
}

/**
 * \brief Evaluate face basis functions for prism cells
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$
 * \param eta2         Reference coordinate \f$\eta_2\f$
 * \param eta3         Reference coordinate \f$\eta_3\f$
 * \param pOrder       Polynomial order
 * \param pos          Face position in reference prism
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives w.r.t. \f$\eta_3\f$
 */
inline void computeFaceBasisAndDerivativesOfPrism
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    const dgFacePositionOnPrism pos,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    switch (pos)
    {
        case dgFacePositionOnPrism::ABE:
            computeInteriorTriBasisAndDerivatives(eta1, eta3, pOrder, basis, dBasis_deta1, dBasis_deta3);
            dBasis_deta2.setSize(basis.size());
            forAll(dBasis_deta2, i)
                dBasis_deta2[i] = 0.0;
            break;
        case dgFacePositionOnPrism::CDF:
            computeInteriorTriBasisAndDerivatives(eta1, eta3, pOrder, basis, dBasis_deta1, dBasis_deta3);
            dBasis_deta2.setSize(basis.size());
            forAll(dBasis_deta2, i)
                dBasis_deta2[i] = 0.0;
            break;
        case dgFacePositionOnPrism::BCEF:
            computeInteriorQuadBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            dBasis_deta1.setSize(basis.size());
            forAll(dBasis_deta1, i)
                dBasis_deta1[i] = 0.0;
            break;
        case dgFacePositionOnPrism::ABCD:
            computeInteriorQuadBasisAndDerivatives(eta1, eta2, pOrder, basis, dBasis_deta1, dBasis_deta2);
            dBasis_deta3.setSize(basis.size());
            forAll(dBasis_deta3, i)
                dBasis_deta3[i] = 0.0;
            break;
        case dgFacePositionOnPrism::ADEF:
            computeInteriorQuadBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            dBasis_deta1.setSize(basis.size());
            forAll(dBasis_deta1, i)
                dBasis_deta1[i] = 0.0;
            break;
        default:
            FatalErrorInFunction
                << "Invalid dgFacePositionOnPrism enum value."
                << abort(FatalError);
    }
}

/**
 * \brief Evaluate face basis functions for pyramid cells
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$
 * \param eta2         Reference coordinate \f$\eta_2\f$
 * \param eta3         Reference coordinate \f$\eta_3\f$
 * \param pOrder       Polynomial order
 * \param pos          Face position in reference pyramid
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives w.r.t. \f$\eta_3\f$
 */
inline void computeFaceBasisAndDerivativesOfPyramid
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    const dgFacePositionOnPyramid pos,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    switch (pos)
    {
        case dgFacePositionOnPyramid::ABCD:
            computeInteriorQuadBasisAndDerivatives(eta1, eta2, pOrder, basis, dBasis_deta1, dBasis_deta2);
            dBasis_deta3.setSize(basis.size());
            forAll(dBasis_deta3, i)
                dBasis_deta3[i] = 0.0;
            break;
        case dgFacePositionOnPyramid::ABE:
            computeInteriorTriBasisAndDerivatives(eta1, eta3, pOrder, basis, dBasis_deta1, dBasis_deta3);
            dBasis_deta2.setSize(basis.size());
            forAll(dBasis_deta2, i)
                dBasis_deta2[i] = 0.0;
            break;
        case dgFacePositionOnPyramid::CDE:
            computeInteriorTriBasisAndDerivatives(eta1, eta3, pOrder, basis, dBasis_deta1, dBasis_deta3);
            dBasis_deta2.setSize(basis.size());
            forAll(dBasis_deta2, i)
                dBasis_deta2[i] = 0.0;
            break;
        case dgFacePositionOnPyramid::ACE:
            computeInteriorTriBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            dBasis_deta1.setSize(basis.size());
            forAll(dBasis_deta1, i)
                dBasis_deta1[i] = 0.0;
            break;
        case dgFacePositionOnPyramid::BDE:
            computeInteriorTriBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            dBasis_deta1.setSize(basis.size());
            forAll(dBasis_deta1, i)
                dBasis_deta1[i] = 0.0;
            break;
        default:
            FatalErrorInFunction
                << "Invalid dgFacePositionOnPyramid enum value."
                << abort(FatalError);
    }
}

/**
 * \brief Evaluate face basis functions for tetrahedral cells
 *
 * \param eta1         Reference coordinate \f$\eta_1\f$
 * \param eta2         Reference coordinate \f$\eta_2\f$
 * \param eta3         Reference coordinate \f$\eta_3\f$
 * \param pOrder       Polynomial order
 * \param pos          Face position in reference tetrahedron
 * \param basis        Basis function values
 * \param dBasis_deta1 Derivatives w.r.t. \f$\eta_1\f$
 * \param dBasis_deta2 Derivatives w.r.t. \f$\eta_2\f$
 * \param dBasis_deta3 Derivatives w.r.t. \f$\eta_3\f$
 */
inline void computeFaceBasisAndDerivativesOfTet
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    const dgFacePositionOnTet pos,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    switch (pos)
    {
        case dgFacePositionOnTet::ABC:
            computeInteriorTriBasisAndDerivatives(eta1, eta2, pOrder, basis, dBasis_deta1, dBasis_deta2);
            break;
        case dgFacePositionOnTet::ABE:
            computeInteriorTriBasisAndDerivatives(eta1, eta3, pOrder, basis, dBasis_deta1, dBasis_deta3);
            break;
        case dgFacePositionOnTet::ACE:
            computeInteriorTriBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            break;
        case dgFacePositionOnTet::BCE:
            computeInteriorTriBasisAndDerivatives(eta2, eta3, pOrder, basis, dBasis_deta2, dBasis_deta3);
            break;
        default:
            FatalErrorInFunction
                << "Invalid dgFacePositionOnTet enum value."
                << abort(FatalError);
    }
}

} // namespace Foam

#endif

// ************************************************************************* //
