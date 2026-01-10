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

#include "dgRefCell.H"
#include "basisFunctions.H"
#include "error.H"
#include "vector.H"
#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dgRefCell::dgRefCell(const label pOrder, const dgCellType type)
:
    pOrder_(pOrder),
    type_(type)
{
    generateCellGaussPointsAndWeights();
    computeBasisAndDerivatives();
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::dgRefCell::getGaussLegendre1D
(
    const label pOrder,
    List<scalar>& eta,
    List<scalar>& w
)
{
    const label nGauss = pOrder + 1;

    eta.setSize(nGauss);
    w.setSize(nGauss);

    switch (nGauss)
    {
        case 1:
            eta[0] = 0.0;
            w [0] = 2.0;
            break;

        case 2:
            eta = { -0.5773502691896257,
                    0.5773502691896257 };
            w  = {  1.0, 1.0 };
            break;

        case 3:
            eta = { -0.7745966692414834,
                    0.0,
                    0.7745966692414834 };
            w  = {  0.5555555555555556,
                    0.8888888888888888,
                    0.5555555555555556 };
            break;

        case 4:
            eta = { -0.8611363115940526,
                   -0.3399810435848563,
                    0.3399810435848563,
                    0.8611363115940526 };
            w  = {  0.3478548451374539,
                    0.6521451548625461,
                    0.6521451548625461,
                    0.3478548451374539 };
            break;

        case 5:
            eta = { -0.9061798459386640,
                   -0.5384693101056831,
                    0.0,
                    0.5384693101056831,
                    0.9061798459386640 };
            w  = {  0.2369268850561891,
                    0.4786286704993665,
                    0.5688888888888889,
                    0.4786286704993665,
                    0.2369268850561891 };
            break;

        case 6:
            eta = { -0.9324695142031521,
                   -0.6612093864662645,
                   -0.2386191860831969,
                    0.2386191860831969,
                    0.6612093864662645,
                    0.9324695142031521 };
            w  = {  0.1713244923791704,
                    0.3607615730481386,
                    0.4679139345726910,
                    0.4679139345726910,
                    0.3607615730481386,
                    0.1713244923791704 };
            break;

        case 7:
            eta = { -0.9491079123427585,
                   -0.7415311855993945,
                   -0.4058451513773972,
                    0.0,
                    0.4058451513773972,
                    0.7415311855993945,
                    0.9491079123427585 };
            w  = {  0.1294849661688697,
                    0.2797053914892766,
                    0.3818300505051189,
                    0.4179591836734694,
                    0.3818300505051189,
                    0.2797053914892766,
                    0.1294849661688697 };
            break;

        default:
            FatalErrorInFunction
                << "Gauss-Legendre not implemented for nGauss = "
                << nGauss << nl
                << "Derived from pOrder = " << pOrder
                << abort(FatalError);
    }
}


void Foam::dgRefCell::getGaussJacobi1D
(
    const label pOrder,
    List<scalar>& eta,
    List<scalar>& w
)
{
    // --------------------------------------------------
    // Fallback conditions
    // --------------------------------------------------
    if (pOrder == 0 || pOrder > 4)
    {
        if (pOrder > 4)
        {
            WarningInFunction
                << "pOrder = " << pOrder
                << " exceeds validated Jacobi range (<=4). "
                << "Falling back to Gauss-Legendre."
                << nl;
        }

        getGaussLegendre1D(pOrder, eta, w);
        return;
    }

    const label nGauss = pOrder + 1;
    eta.setSize(nGauss);
    w.setSize(nGauss);

    switch (nGauss)
    {
        case 2:
            eta = { -0.6546536707079771,
                    0.1546536707079771 };
            w  = {  0.7111111111111111,
                    0.6222222222222222 };
            break;

        case 3:
            eta = { -0.7660444431189780,
                   -0.1736481776669303,
                    0.5 };
            w  = {  0.4444444444444444,
                    0.7777777777777778,
                    0.1111111111111111 };
            break;

        case 4:
            eta = { -0.8302238962785669,
                   -0.4688487934707142,
                    0.0558466708154264,
                    0.7168821974518675 };
            w  = {  0.3038329932404390,
                    0.5909114117464944,
                    0.3719183942939050,
                    0.0666708673851616 };
            break;

        case 5:
            FatalErrorInFunction
                << "Gauss-Jacobi(alpha=2,beta=0) not implemented for nGauss = 5"
                << abort(FatalError);

        default:
            FatalErrorInFunction
                << "Invalid nGauss for Gauss-Jacobi"
                << abort(FatalError);
    }
}


void Foam::dgRefCell::getGaussLobatto1D
(
    const label pOrder,
    List<scalar>& eta,
    List<scalar>& w
)
{
    // --------------------------------------------------
    // Fallback conditions
    // --------------------------------------------------
    if (pOrder == 0 || pOrder > 4)
    {
        if (pOrder > 4)
        {
            WarningInFunction
                << "pOrder = " << pOrder
                << " exceeds validated Lobatto range (<=4). "
                << "Falling back to Gauss-Legendre."
                << nl;
        }

        getGaussLegendre1D(pOrder, eta, w);
        return;
    }

    const label nGauss = pOrder + 1;
    eta.setSize(nGauss);
    w.setSize(nGauss);

    switch (nGauss)
    {
        case 2:
            eta = { -1.0, 1.0 };
            w  = {  1.0, 1.0 };
            break;

        case 3:
            eta = { -1.0, 0.0, 1.0 };
            w  = {  0.3333333333333333,
                    1.3333333333333333,
                    0.3333333333333333 };
            break;

        case 4:
            eta = { -1.0,
                   -0.4472135954999579,
                    0.4472135954999579,
                    1.0 };
            w  = {  0.1666666666666667,
                    0.8333333333333333,
                    0.8333333333333333,
                    0.1666666666666667 };
            break;

        case 5:
            eta = { -1.0,
                   -0.6546536707079771,
                    0.0,
                    0.6546536707079771,
                    1.0 };
            w  = {  0.1,
                    0.5444444444444444,
                    0.7111111111111111,
                    0.5444444444444444,
                    0.1 };
            break;

        default:
            FatalErrorInFunction
                << "Invalid nGauss for Gauss-Lobatto"
                << abort(FatalError);
    }
}


void Foam::dgRefCell::generateCellGaussPointsAndWeights()
{
    // --------------------------------------------------
    // Number of Gauss points per 1D direction
    // Convention: nGauss1D = pOrder + 1
    // --------------------------------------------------
    const label nGauss1D = pOrder_ + 1;

    // --------------------------------------------------
    // 1D Gauss points and weights for each direction
    // --------------------------------------------------
    List<scalar> eta1, w1;
    List<scalar> eta2, w2;
    List<scalar> eta3, w3;

    // --------------------------------------------------
    // Select quadrature rules based on cell type
    // --------------------------------------------------
    switch (type_)
    {
        // ==================================================
        // Hexahedral cell
        // Tensor-product Gauss–Legendre in all directions
        // ==================================================
        case dgCellType::HEX:
        {
            getGaussLegendre1D(pOrder_, eta1, w1);
            getGaussLegendre1D(pOrder_, eta2, w2);
            getGaussLegendre1D(pOrder_, eta3, w3);
            break;
        }

        // ==================================================
        // Tetrahedral cell
        // Collapsed coordinates:
        //   - eta1: Gauss–Jacobi (alpha=2, beta=0)
        //   - eta2, eta3: Gauss–Legendre
        // pOrder = 0 or pOrder > 4 handled internally by fallback
        // ==================================================
        case dgCellType::TET:
        {
            getGaussJacobi1D(pOrder_, eta1, w1);
            getGaussLegendre1D(pOrder_, eta2, w2);
            getGaussLegendre1D(pOrder_, eta3, w3);
            break;
        }

        // ==================================================
        // Prism cell
        // Collapsed direction in eta3:
        //   - eta1, eta2: Gauss–Lobatto
        //   - eta3: Gauss–Radau
        // pOrder = 0 or pOrder > 4 handled internally by fallback
        // ==================================================
        case dgCellType::PRISM:
        {
            getGaussLobatto1D(pOrder_, eta1, w1);
            getGaussLobatto1D(pOrder_, eta2, w2);
            getGaussRadau1D  (pOrder_, eta3, w3);
            break;
        }

        // ==================================================
        // Pyramid cell
        // Same quadrature strategy as prism
        // ==================================================
        case dgCellType::PYRAMID:
        {
            getGaussLobatto1D(pOrder_, eta1, w1);
            getGaussLobatto1D(pOrder_, eta2, w2);
            getGaussRadau1D  (pOrder_, eta3, w3);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unsupported dgCellType"
                << abort(FatalError);
        }
    }

    // --------------------------------------------------
    // Total number of Gauss points (tensor-product)
    // --------------------------------------------------
    nGauss_ = eta1.size()*eta2.size()*eta3.size();

    gaussP_.setSize(nGauss_);
    wGauss_.setSize(nGauss_);

    // --------------------------------------------------
    // Build tensor-product Gauss points and weights
    // Reference coordinates: (eta1, eta2, eta3)
    // --------------------------------------------------
    label idx = 0;

    for (label k = 0; k < eta3.size(); ++k)
    {
        for (label j = 0; j < eta2.size(); ++j)
        {
            for (label i = 0; i < eta1.size(); ++i)
            {
                gaussP_[idx] =
                    vector
                    (
                        eta1[i],
                        eta2[j],
                        eta3[k]
                    );

                wGauss_[idx] =
                    w1[i]*w2[j]*w3[k];

                ++idx;
            }
        }
    }
}


void dgRefCell::computeBasisAndDerivatives()
{
    //const label nBasis = Foam::math::getNumBasis(pOrder_, type_);
    nDof_ = nBasis;

    basis_.setSize(nGauss_);
    dBasis_dEta1_.setSize(nGauss_);
    dBasis_dEta2_.setSize(nGauss_);
    dBasis_dEta3_.setSize(nGauss_);

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        //basis_[gp].setSize(nBasis);
        //dBasis_dEta1_[gp].setSize(nBasis);
        //dBasis_dEta2_[gp].setSize(nBasis);
        //dBasis_dEta3_[gp].setSize(nBasis);

        const scalar eta1 = gaussP_[gp].x();
        const scalar eta2 = gaussP_[gp].y();
        const scalar eta3 = gaussP_[gp].z();

        switch (type_)
        {
            case dgCellType::HEX:
                Foam::math::computeInteriorHexBasisAndDerivatives(
                    eta1, eta2, eta3, pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PRISM:
                Foam::math::computeInteriorPrismBasisAndDerivatives(
                    eta1, eta2, eta3, pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::TET:
                Foam::math::computeInteriorTetBasisAndDerivatives(
                    eta1, eta2, eta3, pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PYRAMID:
                Foam::math::computeInteriorPyramidBasisAndDerivatives(
                    eta1, eta2, eta3, pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            default:
                FatalErrorInFunction
                    << "Unsupported cell type: " << type_ << nl
                    << abort(FatalError);
        }
    }
}

// ************************************************************************* //
