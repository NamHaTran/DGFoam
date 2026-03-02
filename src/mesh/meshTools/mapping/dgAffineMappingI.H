/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#ifndef Foam_dgAffineMappingI_H
#define Foam_dgAffineMappingI_H

#include "dgCellType.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

inline vector mapEtaToXHex(
    const vector& eta, 
    const List<vector>& cellVertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = eta.x();
    const scalar eta2 = eta.y();
    const scalar eta3 = eta.z();

    // ---------------------------------------------------------------------
    // Ordering: A B C D E F G H
    // ---------------------------------------------------------------------

    List<scalar> shapeFunc(8);

    shapeFunc[0] = 0.125 * (1 - eta1) * (1 - eta2) * (1 - eta3); // A
    shapeFunc[1] = 0.125 * (1 + eta1) * (1 - eta2) * (1 - eta3); // B
    shapeFunc[2] = 0.125 * (1 + eta1) * (1 + eta2) * (1 - eta3); // C
    shapeFunc[3] = 0.125 * (1 - eta1) * (1 + eta2) * (1 - eta3); // D
    shapeFunc[4] = 0.125 * (1 - eta1) * (1 - eta2) * (1 + eta3); // E
    shapeFunc[5] = 0.125 * (1 + eta1) * (1 - eta2) * (1 + eta3); // F
    shapeFunc[6] = 0.125 * (1 + eta1) * (1 + eta2) * (1 + eta3); // G
    shapeFunc[7] = 0.125 * (1 - eta1) * (1 + eta2) * (1 + eta3); // H

    vector x(vector::zero);
    forAll(shapeFunc, i)
    {
        x += shapeFunc[i] * cellVertices[i];
    }

    return x;
}

/**
 * \brief Map a physical point to reference coordinates of a trilinear hex
 *
 * \details
 * Solves x = sum_i N_i(eta) Xi using Newton iteration.
 * The mapping is trilinear, therefore non-affine.
 *
 * \param x              Physical point
 * \param cellVertices   List of 8 hex vertices [A,B,C,D,E,F,G,H]
 *
 * \return Reference coordinates (eta1, eta2, eta3)
 */
inline vector mapXToEtaHex
(
    const vector& x,
    const List<vector>& cellVertices
)
{
    const scalar tol = 1e-12;
    const label  maxIter = 20;

    vector eta(0,0,0);   // initial guess (cell center)

    for (label iter = 0; iter < maxIter; ++iter)
    {
        const scalar eta1 = eta.x();
        const scalar eta2 = eta.y();
        const scalar eta3 = eta.z();

        // --------------------------------------------------------------
        // Shape functions
        // --------------------------------------------------------------

        scalar N[8];

        N[0] = 0.125*(1-eta1)*(1-eta2)*(1-eta3);
        N[1] = 0.125*(1+eta1)*(1-eta2)*(1-eta3);
        N[2] = 0.125*(1+eta1)*(1+eta2)*(1-eta3);
        N[3] = 0.125*(1-eta1)*(1+eta2)*(1-eta3);
        N[4] = 0.125*(1-eta1)*(1-eta2)*(1+eta3);
        N[5] = 0.125*(1+eta1)*(1-eta2)*(1+eta3);
        N[6] = 0.125*(1+eta1)*(1+eta2)*(1+eta3);
        N[7] = 0.125*(1-eta1)*(1+eta2)*(1+eta3);

        // --------------------------------------------------------------
        // Compute mapped physical point
        // --------------------------------------------------------------

        vector xEta(vector::zero);

        for (label i = 0; i < 8; ++i)
        {
            xEta += N[i]*cellVertices[i];
        }

        vector F = xEta - x;

        if (mag(F) < tol)
        {
            return eta;
        }

        // --------------------------------------------------------------
        // Compute Jacobian J = dX/deta
        // --------------------------------------------------------------

        tensor J(tensor::zero);

        for (label i = 0; i < 8; ++i)
        {
            const scalar xi = (i==1||i==2||i==5||i==6) ? 1 : -1;
            const scalar yi = ((i>=2 && i<=3) || (i>=6 && i<=7)) ? 1 : -1;
            const scalar zi = (i>=4) ? 1 : -1;

            const scalar dN1 = 0.125*xi*(1+yi*eta2)*(1+zi*eta3);
            const scalar dN2 = 0.125*yi*(1+xi*eta1)*(1+zi*eta3);
            const scalar dN3 = 0.125*zi*(1+xi*eta1)*(1+yi*eta2);

            const vector& Xi = cellVertices[i];

            // Row X
            J.xx() += Xi.x() * dN1;
            J.xy() += Xi.x() * dN2;
            J.xz() += Xi.x() * dN3;

            // Row Y
            J.yx() += Xi.y() * dN1;
            J.yy() += Xi.y() * dN2;
            J.yz() += Xi.y() * dN3;

            // Row Z
            J.zx() += Xi.z() * dN1;
            J.zy() += Xi.z() * dN2;
            J.zz() += Xi.z() * dN3;
        }

        // --------------------------------------------------------------
        // Newton update
        // --------------------------------------------------------------

        if (mag(det(J)) < SMALL)
        {
            FatalErrorInFunction
                << "Singular Jacobian in hex inverse mapping"
                << exit(FatalError);
        }

        vector delta = inv(J) & F;

        eta -= delta;

        if (mag(delta) < tol)
        {
            return eta;
        }
    }

    FatalErrorInFunction
        << "Hex inverse mapping did not converge"
        << exit(FatalError);

    return eta; // never reached
}

inline vector mapEtaToXPrism(
    const vector& eta, 
    const List<vector>& cellVertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = eta.x();
    const scalar eta2 = eta.y();
    const scalar eta3 = eta.z();

    // Vertex correspondence:
    //
    //   cellVertices id -> prism vertex
    //   -------------------------------
    //        0        ->       E       
    //        1        ->       A      
    //        2        ->       B       
    //        3        ->       F       
    //        4        ->       D      
    //        5        ->       C      
    //
    // ---------------------------------------------------------------------

    List<scalar> shapeFunc(6);

    shapeFunc[0] = 0.25 * (1 - eta2) * (1 + eta3); // E
    shapeFunc[1] = 0.125 * (1 - eta1) * (1 - eta2) * (1 - eta3); // A
    shapeFunc[2] = 0.125 * (1 + eta1) * (1 - eta2) * (1 - eta3); // B
    shapeFunc[3] = 0.25 * (1 + eta2) * (1 + eta3); // F
    shapeFunc[4] = 0.125 * (1 - eta1) * (1 + eta2) * (1 - eta3); // D
    shapeFunc[5] = 0.125 * (1 + eta1) * (1 + eta2) * (1 - eta3); // C

    vector x(vector::zero);
    forAll(shapeFunc, i)
    {
        x += shapeFunc[i] * cellVertices[i];
    }

    return x;
}

inline Foam::vector mapXToEtaPyramid
(
    const vector& x,
    const List<vector>& cellVertices
)
{
    const scalar tol     = 1e-12;
    const label  maxIter = 25;

    // Initial guess at reference cell center
    vector eta(0, 0, 0);

    for (label iter = 0; iter < maxIter; ++iter)
    {
        const scalar eta1 = eta.x();
        const scalar eta2 = eta.y();
        const scalar eta3 = eta.z();

        // ==========================================================
        // Shape functions N_i(eta)
        // ==========================================================

        scalar N[5];

        N[0] = 0.25*(1-eta1)*(1-eta2)*(1-eta3); // A
        N[1] = 0.25*(1+eta1)*(1-eta2)*(1-eta3); // B
        N[2] = 0.25*(1+eta1)*(1+eta2)*(1-eta3); // C
        N[3] = 0.25*(1-eta1)*(1+eta2)*(1-eta3); // D
        N[4] = 0.5 *(1+eta3);                   // E (apex)

        // ==========================================================
        // Compute mapped physical point X(eta)
        // ==========================================================

        vector xEta(vector::zero);

        for (label i = 0; i < 5; ++i)
        {
            xEta += N[i] * cellVertices[i];
        }

        vector F = xEta - x;

        // Converged in physical space
        if (mag(F) < tol)
        {
            return eta;
        }

        // ==========================================================
        // Build Jacobian J = dX/deta manually
        //
        // J_ij = sum_k Xk_i * dNk/deta_j
        //
        // Columns:
        //   col 0 -> dX/deta1
        //   col 1 -> dX/deta2
        //   col 2 -> dX/deta3
        // ==========================================================

        tensor J(tensor::zero);

        for (label i = 0; i < 5; ++i)
        {
            scalar dN1, dN2, dN3;

            if (i == 0) // A
            {
                dN1 = -0.25*(1-eta2)*(1-eta3);
                dN2 = -0.25*(1-eta1)*(1-eta3);
                dN3 = -0.25*(1-eta1)*(1-eta2);
            }
            else if (i == 1) // B
            {
                dN1 =  0.25*(1-eta2)*(1-eta3);
                dN2 = -0.25*(1+eta1)*(1-eta3);
                dN3 = -0.25*(1+eta1)*(1-eta2);
            }
            else if (i == 2) // C
            {
                dN1 =  0.25*(1+eta2)*(1-eta3);
                dN2 =  0.25*(1+eta1)*(1-eta3);
                dN3 = -0.25*(1+eta1)*(1+eta2);
            }
            else if (i == 3) // D
            {
                dN1 = -0.25*(1+eta2)*(1-eta3);
                dN2 =  0.25*(1-eta1)*(1-eta3);
                dN3 = -0.25*(1-eta1)*(1+eta2);
            }
            else // E
            {
                dN1 = 0.0;
                dN2 = 0.0;
                dN3 = 0.5;
            }

            const vector& Xi = cellVertices[i];

            // --- Fill J explicitly (no vector*tensor product) ---

            // Row X
            J.xx() += Xi.x() * dN1;
            J.xy() += Xi.x() * dN2;
            J.xz() += Xi.x() * dN3;

            // Row Y
            J.yx() += Xi.y() * dN1;
            J.yy() += Xi.y() * dN2;
            J.yz() += Xi.y() * dN3;

            // Row Z
            J.zx() += Xi.z() * dN1;
            J.zy() += Xi.z() * dN2;
            J.zz() += Xi.z() * dN3;
        }

        // ==========================================================
        // Newton update
        // Solve J * delta = F
        // ==========================================================

        const scalar detJ = det(J);

        if (mag(detJ) < SMALL)
        {
            FatalErrorInFunction
                << "Singular Jacobian in pyramid inverse mapping"
                << exit(FatalError);
        }

        vector delta = inv(J) & F;

        eta -= delta;

        if (mag(delta) < tol)
        {
            return eta;
        }
    }

    FatalErrorInFunction
        << "Pyramid inverse mapping did not converge"
        << exit(FatalError);

    return eta;
}

inline Foam::vector mapXToEtaPrism
(
    const vector& x,
    const List<vector>& cellVertices
)
{
    const scalar tol     = 1e-12;
    const label  maxIter = 20;

    // Initial guess at reference cell center
    vector eta(0, 0, 0);

    for (label iter = 0; iter < maxIter; ++iter)
    {
        const scalar eta1 = eta.x();
        const scalar eta2 = eta.y();
        const scalar eta3 = eta.z();

        // ==========================================================
        // Shape functions
        // ==========================================================

        scalar N[6];

        N[0] = 0.25  * (1 - eta2) * (1 + eta3);                        // E
        N[1] = 0.125 * (1 - eta1) * (1 - eta2) * (1 - eta3);           // A
        N[2] = 0.125 * (1 + eta1) * (1 - eta2) * (1 - eta3);           // B
        N[3] = 0.25  * (1 + eta2) * (1 + eta3);                        // F
        N[4] = 0.125 * (1 - eta1) * (1 + eta2) * (1 - eta3);           // D
        N[5] = 0.125 * (1 + eta1) * (1 + eta2) * (1 - eta3);           // C

        // ==========================================================
        // Compute mapped physical point X(eta)
        // ==========================================================

        vector xEta(vector::zero);

        for (label i = 0; i < 6; ++i)
        {
            xEta += N[i] * cellVertices[i];
        }

        vector F = xEta - x;

        if (mag(F) < tol)
        {
            return eta;
        }

        // ==========================================================
        // Build Jacobian J = dX/deta manually
        //
        // Columns:
        //   col 0 -> dX/deta1
        //   col 1 -> dX/deta2
        //   col 2 -> dX/deta3
        // ==========================================================

        tensor J(tensor::zero);

        for (label i = 0; i < 6; ++i)
        {
            scalar dN1, dN2, dN3;

            if (i == 0) // E
            {
                dN1 = 0.0;
                dN2 = -0.25 * (1 + eta3);
                dN3 =  0.25 * (1 - eta2);
            }
            else if (i == 1) // A
            {
                dN1 = -0.125 * (1 - eta2) * (1 - eta3);
                dN2 = -0.125 * (1 - eta1) * (1 - eta3);
                dN3 = -0.125 * (1 - eta1) * (1 - eta2);
            }
            else if (i == 2) // B
            {
                dN1 =  0.125 * (1 - eta2) * (1 - eta3);
                dN2 = -0.125 * (1 + eta1) * (1 - eta3);
                dN3 = -0.125 * (1 + eta1) * (1 - eta2);
            }
            else if (i == 3) // F
            {
                dN1 = 0.0;
                dN2 =  0.25 * (1 + eta3);
                dN3 =  0.25 * (1 + eta2);
            }
            else if (i == 4) // D
            {
                dN1 = -0.125 * (1 + eta2) * (1 - eta3);
                dN2 =  0.125 * (1 - eta1) * (1 - eta3);
                dN3 = -0.125 * (1 - eta1) * (1 + eta2);
            }
            else // C
            {
                dN1 =  0.125 * (1 + eta2) * (1 - eta3);
                dN2 =  0.125 * (1 + eta1) * (1 - eta3);
                dN3 = -0.125 * (1 + eta1) * (1 + eta2);
            }

            const vector& Xi = cellVertices[i];

            // ---- Fill J explicitly ----

            // Row X
            J.xx() += Xi.x() * dN1;
            J.xy() += Xi.x() * dN2;
            J.xz() += Xi.x() * dN3;

            // Row Y
            J.yx() += Xi.y() * dN1;
            J.yy() += Xi.y() * dN2;
            J.yz() += Xi.y() * dN3;

            // Row Z
            J.zx() += Xi.z() * dN1;
            J.zy() += Xi.z() * dN2;
            J.zz() += Xi.z() * dN3;
        }

        // ==========================================================
        // Newton update
        // Solve J * delta = F
        // ==========================================================

        const scalar detJ = det(J);

        if (mag(detJ) < SMALL)
        {
            FatalErrorInFunction
                << "Singular Jacobian in prism inverse mapping"
                << exit(FatalError);
        }

        vector delta = inv(J) & F;

        eta -= delta;

        if (mag(delta) < tol)
        {
            return eta;
        }
    }

    FatalErrorInFunction
        << "Prism inverse mapping did not converge"
        << exit(FatalError);

    return eta;
}

inline vector mapEtaToXPyramid(
    const vector& eta, 
    const List<vector>& cellVertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = eta.x();
    const scalar eta2 = eta.y();
    const scalar eta3 = eta.z();

    // Vertex correspondence:
    //
    //   cellVertices id -> pyramid vertex
    //   ---------------------------------------------------
    //        0        ->        A         
    //        1        ->        B         
    //        2        ->        C         
    //        3        ->        D         
    //        4        ->        E         
    //
    // ---------------------------------------------------------------------

    List<scalar> shapeFunc(5);

    shapeFunc[0] = 0.25 * (1 - eta1) * (1 - eta2) * (1 - eta3); // A
    shapeFunc[1] = 0.25 * (1 + eta1) * (1 - eta2) * (1 - eta3); // B
    shapeFunc[2] = 0.25 * (1 + eta1) * (1 + eta2) * (1 - eta3); // C
    shapeFunc[3] = 0.25 * (1 - eta1) * (1 + eta2) * (1 - eta3); // D
    shapeFunc[4] = 0.5  * (1 + eta3);                             // E

    vector x(vector::zero);
    forAll(shapeFunc, i)
    {
        x += shapeFunc[i] * cellVertices[i];
    }

    return x;
}

inline vector mapEtaToXTet(
    const vector& eta, 
    const List<vector>& cellVertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = eta.x();
    const scalar eta2 = eta.y();
    const scalar eta3 = eta.z();

    // Vertex correspondence:
    //
    //   cellVertices id -> tetrahedron vertex
    //   ---------------------------------------------------
    //        0        ->        A         
    //        1        ->        B         
    //        2        ->        C         
    //        3        ->        D         
    //
    // ---------------------------------------------------------------------

    List<scalar> shapeFunc(4);

    shapeFunc[0] = 0.125 * (1 - eta1) * (1 - eta2) * (1 - eta3); // A
    shapeFunc[1] = 0.125 * (1 + eta1) * (1 - eta2) * (1 - eta3); // B
    shapeFunc[2] = 0.250 * (1 + eta2) * (1 - eta3); // C
    shapeFunc[3] = 0.5 * (1 + eta3); // D

    vector x(vector::zero);
    forAll(shapeFunc, i)
    {
        x += shapeFunc[i] * cellVertices[i];
    }

    return x;
}

inline vector mapXToEtaTet
(
    const vector& x,
    const List<vector>& cellVertices
)
{
    const scalar eps = SMALL;

    // ---------------------------------------------------------------------
    // Extract tetra vertices
    // ---------------------------------------------------------------------

    const vector& A = cellVertices[0];
    const vector& B = cellVertices[1];
    const vector& C = cellVertices[2];
    const vector& D = cellVertices[3];

    // ---------------------------------------------------------------------
    // Build affine Jacobian matrix M = [B-A, C-A, D-A]
    // Columns of tensor are BA, CA, DA
    // ---------------------------------------------------------------------

    const vector BA = B - A;
    const vector CA = C - A;
    const vector DA = D - A;

    //const tensor M(BA, CA, DA);
    const tensor M
    (
        BA.x(), CA.x(), DA.x(),
        BA.y(), CA.y(), DA.y(),
        BA.z(), CA.z(), DA.z()
    );

    const scalar detM = det(M);

    if (mag(detM) < eps)
    {
        FatalErrorInFunction
            << "Degenerate tetrahedron (determinant ~ 0)"
            << exit(FatalError);
    }

    // ---------------------------------------------------------------------
    // Solve linear system M * lambda = (x - A)
    // ---------------------------------------------------------------------

    const vector rhs = x - A;

    const vector lambda = inv(M) & rhs;

    const scalar lamB = lambda.x();
    const scalar lamC = lambda.y();
    const scalar lamD = lambda.z();
    // const scalar lamA = 1.0 - lamB - lamC - lamD;

    // ---------------------------------------------------------------------
    // Convert barycentric -> collapsed coordinates
    // ---------------------------------------------------------------------

    const scalar eta3 = 2.0*lamD - 1.0;

    const scalar denom2 = 1.0 - lamD;
    if (mag(denom2) < eps)
    {
        FatalErrorInFunction
            << "Singularity near tetra apex"
            << exit(FatalError);
    }

    const scalar eta2 = 2.0*lamC/denom2 - 1.0;

    const scalar denom1 = 1.0 - lamC - lamD;
    if (mag(denom1) < eps)
    {
        FatalErrorInFunction
            << "Singularity near tetra edge"
            << exit(FatalError);
    }

    const scalar eta1 = 2.0*lamB/denom1 - 1.0;

    return vector(eta1, eta2, eta3);
}

inline vector mapEtaToX(
    const vector& eta, 
    const List<vector>& cellVertices,
    const dgCellType cellType
)
{
    switch (cellType)
    {
        case dgCellType::HEX:
            return mapEtaToXHex(eta, cellVertices);
        case dgCellType::PRISM:
            return mapEtaToXPrism(eta, cellVertices);
        case dgCellType::PYRAMID:
            return mapEtaToXPyramid(eta, cellVertices);
        case dgCellType::TET:
            return mapEtaToXTet(eta, cellVertices);
        default:
            FatalErrorInFunction
                << "Unsupported cell shape: " << cellType
                << exit(FatalError);
    }
    return vector::zero; // Should never reach here
}

inline vector mapXToEta(
    const vector& x, 
    const List<vector>& cellVertices,
    const dgCellType cellType
)
{
    switch (cellType)
    {
        case dgCellType::HEX:
            return mapXToEtaHex(x, cellVertices);
        case dgCellType::PRISM:
            return mapXToEtaPrism(x, cellVertices);
        case dgCellType::PYRAMID:
            return mapXToEtaPyramid(x, cellVertices);
        case dgCellType::TET:
            return mapXToEtaTet(x, cellVertices);
        default:
            FatalErrorInFunction
                << "Unsupported cell shape: " << cellType
                << exit(FatalError);
    }
    return vector::zero; // Should never reach here
}

inline vector mapEtaToXiPrism(
    const vector& eta 
)
{
    // Extract prism coordinates
    const scalar eta1 = eta.x();
    const scalar eta2 = eta.y();
    const scalar eta3 = eta.z();

    // xi3 equals eta3
    const scalar xi3 = eta3;

    // Scale factor
    const scalar scale = 1.0 - eta3;

    const scalar xi1 = 0.5*(1.0 + eta1)*scale - 1.0;
    const scalar xi2 = eta2;

    return vector(xi1, xi2, xi3);
}

inline vector mapXiToEtaPrism(
    const vector& xi 
)
{
    // Extract hex coordinates
    const scalar xi1 = xi.x();
    const scalar xi2 = xi.y();
    const scalar xi3 = xi.z();

    // Denominator (collapse direction)
    const scalar denom = 1.0 - xi3;

    // NOTE:
    // xi3 = 1 causes singularity (collapsed triangle apex)
    // Ensure quadrature does not include xi3 = 1

    const scalar eta1 = 2.0*(1.0 + xi1)/denom - 1.0;
    const scalar eta2 = xi2;
    const scalar eta3 = xi3;

    return vector(eta1, eta2, eta3);
}

inline vector mapEtaToXiPyramid(
    const vector& eta 
)
{
    // Extract eta coordinates
    const scalar eta1 = eta.x();
    const scalar eta2 = eta.y();
    const scalar eta3 = eta.z();

    // xi3 equals eta3 directly
    const scalar xi3 = eta3;

    // Scaling factor (1 - eta3)
    const scalar scale = 1.0 - eta3;

    const scalar xi1 = 0.5*(1.0 + eta1)*scale - 1.0;
    const scalar xi2 = 0.5*(1.0 + eta2)*scale - 1.0;

    return vector(xi1, xi2, xi3);
}

inline vector mapXiToEtaPyramid(
    const vector& xi 
)
{
    // Extract xi coordinates
    const scalar xi1 = xi.x();
    const scalar xi2 = xi.y();
    const scalar xi3 = xi.z();

    // Denominator term (1 - xi3)
    const scalar denom = 1.0 - xi3;

    // NOTE:
    // xi3 = 1 is singular (collapsed top of prism)
    // User must ensure xi3 != 1 in quadrature

    const scalar eta1 = 2.0*(1.0 + xi1)/denom - 1.0;
    const scalar eta2 = 2.0*(1.0 + xi2)/denom - 1.0;
    const scalar eta3 = xi3;

    return vector(eta1, eta2, eta3);
}

// Map from reference hex (xi) to reference tetra (eta)
inline vector mapXiToEtaTetra(const vector& xi)
{
    // Extract hex coordinates
    const scalar xi1 = xi.x();
    const scalar xi2 = xi.y();
    const scalar xi3 = xi.z();

    // Denominators of the Duffy transform
    const scalar denom1 = -xi2 - xi3;  // collapses (xi2, xi3) to an edge
    const scalar denom2 =  1.0 - xi3;  // collapses xi3 to a face

    // NOTE:
    // denom1 = 0 or denom2 = 0 makes the mapping singular.
    // Ensure quadrature does not include these locations.

    const scalar eta1 = 2.0*(1.0 + xi1)/denom1 - 1.0;
    const scalar eta2 = 2.0*(1.0 + xi2)/denom2 - 1.0;
    const scalar eta3 = xi3;

    return vector(eta1, eta2, eta3);
}


// Map from reference tetra (eta) back to reference hex (xi)
inline vector mapEtaToXiTetra(const vector& eta)
{
    // Extract tetra coordinates
    const scalar eta1 = eta.x();
    const scalar eta2 = eta.y();
    const scalar eta3 = eta.z();

    // Direct mapping for the third coordinate
    const scalar xi3 = eta3;

    // Recover xi2 from eta2 equation
    const scalar xi2 = 0.5*(1.0 + eta2)*(1.0 - eta3) - 1.0;

    // Recover xi1 (expanded inverse form; avoids using xi2 in denom explicitly)
    const scalar xi1 = 0.25*(1.0 + eta1)*(1.0 - eta2)*(1.0 - eta3) - 1.0;

    return vector(xi1, xi2, xi3);
}

inline Foam::vector mapToPhysicalTriangle
(
    const vector2D& rs,
    const List<vector>& vertices
)
{
    const scalar r = rs.x();
    const scalar s = rs.y();

    return (1.0 - r - s)*vertices[0]
         + r*vertices[1]
         + s*vertices[2];
}

inline Foam::vector mapToPhysicalQuad
(
    const vector2D& eta,
    const List<vector>& vertices
)
{
    const scalar eta1  = eta.x();
    const scalar eta2 = eta.y();

    const scalar N0 = 0.25*(1 - eta1)*(1 - eta2);
    const scalar N1 = 0.25*(1 + eta1)*(1 - eta2);
    const scalar N2 = 0.25*(1 + eta1)*(1 + eta2);
    const scalar N3 = 0.25*(1 - eta1)*(1 + eta2);

    return N0*vertices[0]
         + N1*vertices[1]
         + N2*vertices[2]
         + N3*vertices[3];
}

} // End namespace Foam

#endif

// ************************************************************************* //
