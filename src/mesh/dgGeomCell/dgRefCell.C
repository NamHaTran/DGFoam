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

void dgRefCell::generateCellGaussPointsAndWeights()
{
    List<scalar> eta1D, w1D;

    switch (pOrder_)
    {
        case 0:
            eta1D = { 0.0 };
            w1D   = { 2.0 };
            break;

        case 1:
            eta1D = { -0.5773502692, 0.5773502692 };
            w1D   = { 1.0, 1.0 };
            break;

        case 2:
            eta1D = { -0.7745966692, 0.0, 0.7745966692 };
            w1D   = { 0.5555555556, 0.8888888889, 0.5555555556 };
            break;

        case 3:
            eta1D = { -0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116 };
            w1D   = { 0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451 };
            break;

        default:
            FatalErrorInFunction
                << "Gauss rule for pOrder = " << pOrder_ << " not implemented\n"
                << abort(FatalError);
    }

    const label n1D = eta1D.size();
    nGauss_ = n1D * n1D * n1D;

    gaussP_.setSize(nGauss_);
    wGauss_.setSize(nGauss_);

    label idx = 0;
    for (label k = 0; k < n1D; ++k)
    {
        for (label j = 0; j < n1D; ++j)
        {
            for (label i = 0; i < n1D; ++i)
            {
                gaussP_[idx] = vector(eta1D[i], eta1D[j], eta1D[k]); // (eta1, eta2, eta3)
                wGauss_[idx] = w1D[i] * w1D[j] * w1D[k];
                ++idx;
            }
        }
    }
}


void dgRefCell::computeBasisAndDerivatives()
{
    const label nBasis = Foam::math::getNumBasis(pOrder_, type_);
    nDof_ = nBasis;

    basis_.setSize(nGauss_);
    dBasis_dEta1_.setSize(nGauss_);
    dBasis_dEta2_.setSize(nGauss_);
    dBasis_dEta3_.setSize(nGauss_);

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        basis_[gp].setSize(nBasis);
        dBasis_dEta1_[gp].setSize(nBasis);
        dBasis_dEta2_[gp].setSize(nBasis);
        dBasis_dEta3_[gp].setSize(nBasis);

        const scalar eta1 = gaussP_[gp].x();
        const scalar eta2 = gaussP_[gp].y();
        const scalar eta3 = gaussP_[gp].z();

        switch (type_)
        {
            case dgCellType::HEX:
                Foam::math::computeHexBasisAndDerivatives(
                    eta1, eta2, eta3, pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PRISM:
                Foam::math::computePrismBasisAndDerivatives(
                    eta1, eta2, eta3, pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::TET:
                Foam::math::computeTetBasisAndDerivatives(
                    eta1, eta2, eta3, pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PYRAMID:
                Foam::math::computePyramidBasisAndDerivatives(
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
