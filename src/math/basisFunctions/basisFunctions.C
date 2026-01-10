/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "basisFunctions.H"
#include "error.H"
#include <cmath>
#include <boost/math/special_functions/jacobi.hpp>
#include <boost/math/tools/config.hpp>
#include "scalar.H"
#include "label.H"
#include "List.H"
#include "dgFacePosition.H"

namespace Foam
{
namespace math
{

label getNumBasis(const label pOrder, const dgCellType type)
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

/*
void computeBasisAndDerivatives
(
    const scalar eta1,
    const scalar eta2,
    const scalar eta3,
    const label pOrder,
    const dgCellType type,
    List<scalar>& basis,
    List<scalar>& dBasis_deta1,
    List<scalar>& dBasis_deta2,
    List<scalar>& dBasis_deta3
)
{
    switch (type)
    {
        case dgCellType::TET:
            // computeTetBasis(eta1, eta2, eta3, pOrder, basis);
            break;

        case dgCellType::HEX:
            // computeHexBasis(eta1, eta2, eta3, pOrder, basis);
            break;

        case dgCellType::PRISM:
            // computePrismBasis(eta1, eta2, eta3, pOrder, basis);
            break;

        case dgCellType::PYRAMID:
            // computePyramidBasis(eta1, eta2, eta3, pOrder, basis);
            break;
        case dgCellType::INVALID:
        default:
            FatalErrorInFunction
                << "Unsupported or invalid cell type for basis function: "
                << static_cast<int>(type) << nl
                << exit(FatalError);
    }
}
*/

void computeInteriorHexBasisAndDerivatives
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

    const label nBasis = getNumBasis(pOrder, dgCellType::HEX);
    basis.setSize(nBasis);
    dBasis_deta1.setSize(nBasis);
    dBasis_deta2.setSize(nBasis);
    dBasis_deta3.setSize(nBasis);

    label idx = 0;

    label polyDegree = pOrder+1;

    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp   = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp  = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);

        for (label q = 1; q < polyDegree; ++q)
        {
            scalar Pq   = jacobi(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);
            scalar dPq  = jacobi_prime(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);

            for (label r = 1; r < polyDegree; ++r)
            {
                scalar Pr   = jacobi(static_cast<unsigned>(r-1), 1.0, 1.0, eta3);
                scalar dPr  = jacobi_prime(static_cast<unsigned>(r-1), 1.0, 1.0, eta3);

                basis[idx]        = Pp * Pq * Pr;
                dBasis_deta1[idx] = dPp * Pq * Pr;
                dBasis_deta2[idx] = Pp  * dPq * Pr;
                dBasis_deta3[idx] = Pp  * Pq  * dPr;

                ++idx;
            }
        }
    }
}


void computeInteriorPrismBasisAndDerivatives
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


void computeInteriorPyramidBasisAndDerivatives
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


void computeInteriorTetBasisAndDerivatives
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

void computeInteriorTriBasisAndDerivatives
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
    // using boost::math::jacobi_prime;

    basis.setSize((pOrder-1)*pOrder/2);

    label polyDegree = pOrder+1;
    idx = 0;

    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        
        const scalar oneMinusEta2 = 1.0 - eta2;
        const scalar oneMinusEta2PowP = std::pow(oneMinusEta2, p);
        const scalar dOneMinusEta2PowP = (p > 0)
            ? -p * std::pow(oneMinusEta2, p - 1)
            : 0.0;
        
        const double alpha = 2 * q + 1;

        for (label q = 1; q < polyDegree; ++q)
        {
            if (p + q < polyDegree)
            {
                scalar Pq  = jacobi(static_cast<unsigned>(q-1), alpha, 1.0, eta2);
                // Basis value
                basis[idx] = Pp * oneMinusEta2PowP * Pq;
                dBasis_deta1[idx] = dPp * oneMinusEta2PowP * Pq;
                dBasis_deta2[idx] = Pp *
                    (dOneMinusEta2PowP * Pq + oneMinusEta2PowP * dPq);
                
                ++idx;
            }
        }
    }
}

void computeInteriorQuadBasisAndDerivatives
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
    // using boost::math::jacobi_prime;

    basis.setSize(pOrder*pOrder/2);

    label polyDegree = pOrder+1;
    idx = 0;
    
    for (label p = 1; p < polyDegree; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);
        scalar dPp    = jacobi_prime(static_cast<unsigned>(p-1), 1.0, 1.0, eta1);

        for (label q = 1; q < polyDegree; ++q)
        {
            scalar Pq  = jacobi(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);
            scalar dPq    = jacobi_prime(static_cast<unsigned>(q-1), 1.0, 1.0, eta2);
            // Basis value
            basis[idx] = Pp * Pq;
            dBasis_deta1[idx] = dPp * Pq;
            dBasis_deta2[idx] = Pp * dPq;

            ++idx;
        }
    }
}

void computeFaceBasisOfHex
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

void computeFaceBasisAndDerivativesOfPrism
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

void computeFaceBasisAndDerivativesOfPyramid
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
    ABCD = 0,
    ACE,
    CDE,
    BDE,
    ABE,

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

void computeFaceBasisAndDerivativesOfTet
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
            computeInteriorTriBasisAndDerivatives(eta1, eta2, pOrder, basis);
            break;
        case dgFacePositionOnTet::ABE:
            computeInteriorTriBasisAndDerivatives(eta1, eta3, pOrder, basis);
            break;
        case dgFacePositionOnTet::ACE:
            computeInteriorTriBasisAndDerivatives(eta2, eta3, pOrder, basis);
            break;
        case dgFacePositionOnTet::BCE:
            computeInteriorTriBasisAndDerivatives(eta2, eta3, pOrder, basis);
            break;
        default:
            FatalErrorInFunction
                << "Invalid dgFacePositionOnTet enum value."
                << abort(FatalError);
    }
}

} // namespace math
} // namespace Foam
// ************************************************************************* //
