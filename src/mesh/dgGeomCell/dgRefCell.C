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
#include "Jacobian.H"
#include "error.H"
#include "vector.H"
#include <cmath>
#include "GaussQuadrature.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dgRefCell::dgRefCell(const label pOrder, const dgCellType type)
:
    pOrder_(pOrder),
    type_(type)
{
    generateCellGaussPointsAndWeights();
    computeBasisAndDerivatives();
    constructMassMatrix();
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::dgRefCell::generateCellGaussPointsAndWeights()
{
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
            Foam::getGaussLobatto1D(pOrder_, eta1, w1);
            Foam::getGaussLobatto1D(pOrder_, eta2, w2);
            Foam::getGaussLobatto1D(pOrder_, eta3, w3);
            break;
        }

        // ==================================================
        // Tetrahedral cell
        // Collapsed coordinates:
        //   - eta1: Gauss–Jacobi (alpha=2, beta=0)
        //   - eta2, eta3: Gauss–Legendre
        // pOrder = 0 handled internally by fallback
        // ==================================================
        case dgCellType::TET:
        {
            Foam::getGaussLobatto1D(pOrder_, eta1, w1);
            Foam::getGaussRadau1D(pOrder_, true, eta2, w2);
            Foam::getGaussRadau1D(pOrder_, true, eta3, w3);
            break;
        }

        // ==================================================
        // Prism cell
        // Collapsed direction in eta3:
        //   - eta1, eta2: Gauss–Lobatto
        //   - eta3: Gauss–Radau
        // pOrder = 0 handled internally by fallback
        // ==================================================
        case dgCellType::PRISM:
        {
            Foam::getGaussLobatto1D(pOrder_, eta1, w1);
            Foam::getGaussLobatto1D(pOrder_, eta2, w2);
            Foam::getGaussRadau1D  (pOrder_, true, eta3, w3);
            break;
        }

        // ==================================================
        // Pyramid cell
        // Same quadrature strategy as prism
        // ==================================================
        case dgCellType::PYRAMID:
        {
            Foam::getGaussLobatto1D(pOrder_, eta1, w1);
            Foam::getGaussLobatto1D(pOrder_, eta2, w2);
            Foam::getGaussRadau1D  (pOrder_, true, eta3, w3);
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
    const label nBasis = Foam::getNumBasis(pOrder_, type_);
    nDof_ = nBasis;

    basis_.setSize(nGauss_);
    dBasis_dEta1_.setSize(nGauss_);
    dBasis_dEta2_.setSize(nGauss_);
    dBasis_dEta3_.setSize(nGauss_);

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        switch (type_)
        {
            case dgCellType::HEX:
                Foam::computeHexBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PRISM:
                Foam::computePrismBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::TET:
                Foam::computeTetBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PYRAMID:
                Foam::computePyramidBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
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

void Foam::dgRefCell::constructMassMatrix()
{
    // Resize and initialize mass matrix
    massMatrix_.resize(nDof_);

    // Resize mass matrix diagonal
    massMatrixDiag_.setSize(nDof_);

    // Assemble mass matrix directly at Gauss points
    for (label i = 0; i < nDof_; ++i)
    {
        for (label j = 0; j < nDof_; ++j)
        {
            scalar mij = 0.0;

            for (label gp = 0; gp < nGauss_; ++gp)
            {
                scalar JDuffy = 0.0;

                switch (type_)
                {
                    case dgCellType::HEX:
                        JDuffy = 1.0; // No Duffy transformation for hex
                        break;

                    case dgCellType::PRISM:
                        JDuffy = Foam::referenceJacobian::prismRefToHexRef(gaussP_[gp]);
                        break;

                    case dgCellType::PYRAMID:
                        JDuffy = Foam::referenceJacobian::pyramidRefToHexRef(gaussP_[gp]);
                        break;

                    case dgCellType::TET:
                        JDuffy = Foam::referenceJacobian::tetRefToHexRef(gaussP_[gp]);
                        break;

                    default:
                        FatalErrorInFunction
                            << "Unsupported cell type for mass matrix calculation: " << type_ << nl
                            << abort(FatalError);
                }

                mij +=
                    basis_[gp][i]
                  * basis_[gp][j]
                  * JDuffy
                  * wGauss_[gp];
                
                // Store diagonal entries separately for efficient mass matrix inversion
                if (i == j)
                {
                    massMatrixDiag_[i] = mij;
                }
            }

            massMatrix_[i][j] = mij;
        }
    }
}

// ************************************************************************* //
