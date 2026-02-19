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

#include "dgGeomFace.H"
#include "pointField.H"
#include "polyMesh.H"

#include "dgRefFace.H"
#include "basisFunctions.H"
#include "error.H"
#include "vector.H"
#include <cmath>

// * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

// (Không có ở đây, thêm sau nếu cần)

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dgRefFace::dgRefFace(const label pOrder)
:
    pOrder_(pOrder)
{
    generateFaceGaussPointsAndWeights();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dgRefFace::generateFaceGaussPointsAndWeights()
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

        case 4:
            eta1D = { -0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459 };
            w1D   = { 0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851 };
            break;

        default:
            FatalErrorInFunction
                << "Gauss rule for pOrder = " << pOrder_ << " not implemented\n"
                << abort(FatalError);
    }

    label n1D = eta1D.size();
    nGauss_ = n1D * n1D;

    gaussP_ABCD_.setSize(nGauss_);
    gaussP_EFGH_.setSize(nGauss_);
    gaussP_ABEF_.setSize(nGauss_);
    gaussP_CDGH_.setSize(nGauss_);
    gaussP_BDFH_.setSize(nGauss_);
    gaussP_ACEG_.setSize(nGauss_);
    wGauss_.setSize(nGauss_);

    label idx = 0;
    for (label j = 0; j < n1D; ++j)
    {
        for (label i = 0; i < n1D; ++i)
        {
            scalar eta1 = eta1D[i];
            scalar eta2 = eta1D[j];

            gaussP_ABCD_[idx] = vector(eta1, eta2, -1.0);
            gaussP_EFGH_[idx] = vector(eta1, eta2,  1.0);
            gaussP_ABEF_[idx] = vector(eta1, -1.0, eta2);
            gaussP_CDGH_[idx] = vector(eta1,  1.0, eta2);
            gaussP_ACEG_[idx] = vector(-1.0, eta1, eta2);
            gaussP_BDFH_[idx] = vector( 1.0, eta1, eta2);

            wGauss_[idx] = w1D[i] * w1D[j];
            ++idx;
        }
    }
}


const Foam::List<Foam::vector>& Foam::dgRefFace::gaussPoints(const dgFacePosition pos) const
{
    switch (pos)
    {
        case dgFacePosition::ABCD:
            return gaussP_ABCD_;
        case dgFacePosition::EFGH:
            return gaussP_EFGH_;
        case dgFacePosition::ABEF:
            return gaussP_ABEF_;
        case dgFacePosition::CDGH:
            return gaussP_CDGH_;
        case dgFacePosition::ACEG:
            return gaussP_ACEG_;
        case dgFacePosition::BDFH:
            return gaussP_BDFH_;
        default:
            FatalErrorInFunction
                << "Invalid dgFacePosition enum value."
                << abort(FatalError);
    }

    // Dummy
    return gaussP_ABCD_; // Should never reach here
}

Foam::basisData Foam::dgRefFace::computeBasisAndDerivatives
(
    const dgCellType cellType, // Cell type for basis functions
    const dgFacePosition pos // Face position for which to compute basis functions
)
{
    // const label nBasis = Foam::math::getNumBasis(pOrder_, cellType);

    Foam::basisData basisData;
    basisData.basis.setSize(nGauss_);
    basisData.dBasis_dEta1.setSize(nGauss_);
    basisData.dBasis_dEta2.setSize(nGauss_);
    basisData.dBasis_dEta3.setSize(nGauss_);

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        vector etaPt;

        // Get Gauss point coordinates based on face position
        // and compute basis functions and derivatives
        switch (pos)
        {
            case dgFacePosition::ABCD:
                etaPt = gaussP_ABCD_[gp];
                break;
            case dgFacePosition::EFGH:
                etaPt = gaussP_EFGH_[gp];
                break;
            case dgFacePosition::ABEF:
                etaPt = gaussP_ABEF_[gp];
                break;
            case dgFacePosition::CDGH:
                etaPt = gaussP_CDGH_[gp];
                break;
            case dgFacePosition::ACEG:
                etaPt = gaussP_ACEG_[gp];
                break;
            case dgFacePosition::BDFH:
                etaPt = gaussP_BDFH_[gp];
                break;
            default:
                FatalErrorInFunction
                    << "Invalid dgFacePosition enum value."
                    << abort(FatalError);
        }
        
        switch (cellType)
        {
            case dgCellType::HEX:
                Foam::computeHexBasisAndDerivatives(
                    etaPt, pOrder_,
                    //pos,
                    basisData.basis[gp],
                    basisData.dBasis_dEta1[gp],
                    basisData.dBasis_dEta2[gp],
                    basisData.dBasis_dEta3[gp]
                );
                break;

            case dgCellType::PRISM:
                Foam::computePrismBasisAndDerivatives(
                    etaPt, pOrder_,
                    //mapFacePositionToPrism(pos),
                    basisData.basis[gp],
                    basisData.dBasis_dEta1[gp],
                    basisData.dBasis_dEta2[gp],
                    basisData.dBasis_dEta3[gp]
                );
                break;

            case dgCellType::TET:
                Foam::computeTetBasisAndDerivatives(
                    etaPt, pOrder_,
                    //mapFacePositionToTet(pos),
                    basisData.basis[gp],
                    basisData.dBasis_dEta1[gp],
                    basisData.dBasis_dEta2[gp],
                    basisData.dBasis_dEta3[gp]
                );
                break;

            case dgCellType::PYRAMID:
                Foam::computePyramidBasisAndDerivatives(
                    etaPt, pOrder_,
                    //mapFacePositionToPyramid(pos),
                    basisData.basis[gp],
                    basisData.dBasis_dEta1[gp],
                    basisData.dBasis_dEta2[gp],
                    basisData.dBasis_dEta3[gp]
                );
                break;

            default:
                FatalErrorInFunction
                    << "Unsupported cell type: " << cellType << nl
                    << abort(FatalError);
        }
    }

    return basisData;
}
// ************************************************************************* //

