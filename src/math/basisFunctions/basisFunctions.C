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
            for (label p = 0; p <= pOrder; ++p)
                for (label q = 0; q <= pOrder; ++q)
                    for (label r = 0; r <= pOrder; ++r)
                        ++count;
            return count;
        }

        case dgCellType::PRISM:
        {
            for (label p = 0; p <= pOrder; ++p)
                for (label q = 0; q <= pOrder; ++q)
                    for (label r = 0; r <= pOrder; ++r)
                        if (q + r <= pOrder)
                            ++count;
            return count;
        }

        case dgCellType::PYRAMID:
        {
            for (label p = 0; p <= pOrder; ++p)
                for (label q = 0; q <= pOrder; ++q)
                    for (label r = 0; r <= pOrder; ++r)
                        if (p + q + r <= pOrder)
                            ++count;
            return count;
        }

        case dgCellType::TET:
        {
            for (label p = 0; p <= pOrder; ++p)
                for (label q = 0; q <= pOrder - p; ++q)
                    for (label r = 0; r <= pOrder - p - q; ++r)
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
            return -1; // Just to satisfy compiler
        }
    }
}


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


void computeHexBasisAndDerivatives
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

    for (label p = 0; p <= pOrder; ++p)
    {
        scalar Pp   = jacobi(static_cast<unsigned>(p), 0.0, 0.0, eta1);
        scalar dPp  = jacobi_prime(static_cast<unsigned>(p), 0.0, 0.0, eta1);

        for (label q = 0; q <= pOrder; ++q)
        {
            scalar Pq   = jacobi(static_cast<unsigned>(q), 0.0, 0.0, eta2);
            scalar dPq  = jacobi_prime(static_cast<unsigned>(q), 0.0, 0.0, eta2);

            for (label r = 0; r <= pOrder; ++r)
            {
                scalar Pr   = jacobi(static_cast<unsigned>(r), 0.0, 0.0, eta3);
                scalar dPr  = jacobi_prime(static_cast<unsigned>(r), 0.0, 0.0, eta3);

                basis[idx]        = Pp * Pq * Pr;
                dBasis_deta1[idx] = dPp * Pq * Pr;
                dBasis_deta2[idx] = Pp  * dPq * Pr;
                dBasis_deta3[idx] = Pp  * Pq  * dPr;

                ++idx;
            }
        }
    }
}


void computePrismBasisAndDerivatives
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

    for (label p = 0; p <= pOrder; ++p)
    {
        scalar Pp     = jacobi(static_cast<unsigned>(p), 0.0, 0.0, eta1);
        scalar dPp    = jacobi_prime(static_cast<unsigned>(p), 0.0, 0.0, eta1);

        for (label q = 0; q <= pOrder; ++q)
        {
            scalar Pq     = jacobi(static_cast<unsigned>(q), 0.0, 0.0, eta2);
            scalar dPq    = jacobi_prime(static_cast<unsigned>(q), 0.0, 0.0, eta2);

            const scalar oneMinusEta3 = 1.0 - eta3;
            const scalar oneMinusEta3PowP = std::pow(oneMinusEta3, p);
            const scalar dOneMinusEta3PowP = (p > 0)
                ? -p * std::pow(oneMinusEta3, p - 1)
                : 0.0;

            const double alpha = 2 * q + 1;

            for (label r = 0; r <= pOrder; ++r)
            {
                if (q + r <= pOrder)
                {
                    scalar Pr    = jacobi(static_cast<unsigned>(r), alpha, 0.0, eta3);
                    scalar dPr   = jacobi_prime(static_cast<unsigned>(r), alpha, 0.0, eta3);

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


void computePyramidBasisAndDerivatives
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

    for (label p = 0; p <= pOrder; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p), 0.0, 0.0, eta1);
        scalar dPp = jacobi_prime(static_cast<unsigned>(p), 0.0, 0.0, eta1);

        for (label q = 0; q <= pOrder; ++q)
        {
            scalar Pq  = jacobi(static_cast<unsigned>(q), 0.0, 0.0, eta2);
            scalar dPq = jacobi_prime(static_cast<unsigned>(q), 0.0, 0.0, eta2);

            double exponent = p + q;
            scalar oneMinusEta3 = 1.0 - eta3;
            scalar powTerm = std::pow(oneMinusEta3, exponent);
            scalar dPowTerm = (exponent > 0)
                ? -exponent * std::pow(oneMinusEta3, exponent - 1)
                : 0.0;

            for (label r = 0; r <= pOrder; ++r)
            {
                if (p + q + r <= pOrder)
                {
                    double alpha = 2 * p + 2 * q + 1;
                    scalar Pr   = jacobi(static_cast<unsigned>(r), alpha, 0.0, eta3);
                    scalar dPr  = jacobi_prime(static_cast<unsigned>(r), alpha, 0.0, eta3);

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


void computeTetBasisAndDerivatives
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

    for (label p = 0; p <= pOrder; ++p)
    {
        scalar Pp  = jacobi(static_cast<unsigned>(p), 0.0, 0.0, eta1);
        scalar dPp = jacobi_prime(static_cast<unsigned>(p), 0.0, 0.0, eta1);

        scalar oneMinusEta2 = 1.0 - eta2;
        scalar oneMinusEta2PowP = std::pow(oneMinusEta2, p);
        scalar dOneMinusEta2PowP = (p > 0) ? -p * std::pow(oneMinusEta2, p - 1) : 0.0;

        for (label q = 0; q <= pOrder - p; ++q)
        {
            double alpha_q = 2 * p + 1;
            scalar Pq     = jacobi(static_cast<unsigned>(q), alpha_q, 0.0, eta2);
            scalar dPq    = jacobi_prime(static_cast<unsigned>(q), alpha_q, 0.0, eta2);

            scalar oneMinusEta3 = 1.0 - eta3;
            double pqSum = p + q;
            scalar oneMinusEta3PowPQ = std::pow(oneMinusEta3, pqSum);
            scalar dOneMinusEta3PowPQ = (pqSum > 0)
                ? -pqSum * std::pow(oneMinusEta3, pqSum - 1)
                : 0.0;

            for (label r = 0; r <= pOrder - p - q; ++r)
            {
                double alpha_r = 2 * p + 2 * q + 1;
                scalar Pr  = jacobi(static_cast<unsigned>(r), alpha_r, 0.0, eta3);
                scalar dPr = jacobi_prime(static_cast<unsigned>(r), alpha_r, 0.0, eta3);

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


} // namespace math
} // namespace Foam
// ************************************************************************* //
