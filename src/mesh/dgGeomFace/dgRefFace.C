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
#include "GaussQuadrature.H"

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

    // Use 1D Gauss-Legendre quadrature for face integration
    Foam::getGaussLegendre1D(pOrder_, eta1D, w1D);

    label n1D = eta1D.size();
    nGauss_ = n1D * n1D;

    gaussP_2D_.setSize(nGauss_);
    wGauss_.setSize(nGauss_);

    label idx = 0;
    for (label j = 0; j < n1D; ++j)
    {
        for (label i = 0; i < n1D; ++i)
        {
            scalar eta1 = eta1D[i];
            scalar eta2 = eta1D[j];

            gaussP_2D_[idx] = vector2D(eta1, eta2);

            wGauss_[idx] = w1D[i] * w1D[j];
            ++idx;
        }
    }
}

Foam::basisData Foam::dgRefFace::computeBasisAndDerivatives
(
    const List<vector>& etaPt,// List ofGauss point in reference coordinates
    const dgCellType cellType // Cell type for basis functions
)
{
    Foam::basisData basisData;
    basisData.basis.setSize(nGauss_);
    basisData.dBasis_dEta1.setSize(nGauss_);
    basisData.dBasis_dEta2.setSize(nGauss_);
    basisData.dBasis_dEta3.setSize(nGauss_);

    for (label gp = 0; gp < nGauss_; ++gp)
    {   
        switch (cellType)
        {
            case dgCellType::HEX:
                Foam::computeHexBasisAndDerivatives(
                    etaPt[gp], pOrder_,
                    //pos,
                    basisData.basis[gp],
                    basisData.dBasis_dEta1[gp],
                    basisData.dBasis_dEta2[gp],
                    basisData.dBasis_dEta3[gp]
                );
                break;

            case dgCellType::PRISM:
                Foam::computePrismBasisAndDerivatives(
                    etaPt[gp], pOrder_,
                    //mapFacePositionToPrism(pos),
                    basisData.basis[gp],
                    basisData.dBasis_dEta1[gp],
                    basisData.dBasis_dEta2[gp],
                    basisData.dBasis_dEta3[gp]
                );
                break;

            case dgCellType::TET:
                Foam::computeTetBasisAndDerivatives(
                    etaPt[gp], pOrder_,
                    //mapFacePositionToTet(pos),
                    basisData.basis[gp],
                    basisData.dBasis_dEta1[gp],
                    basisData.dBasis_dEta2[gp],
                    basisData.dBasis_dEta3[gp]
                );
                break;

            case dgCellType::PYRAMID:
                Foam::computePyramidBasisAndDerivatives(
                    etaPt[gp], pOrder_,
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

Foam::scalar Foam::dgRefFace::calcDuffyJacobian
(
    const dgCellType cellType,
    const dgFacePosition facePos,
    const vector& gaussPt
)
{
    Foam::scalar JDuffy = 1;
    switch (cellType)
    {
        case dgCellType::HEX:
            JDuffy = 1.0/4.0; // No Duffy transformation needed for hexahedron
            break;

        case dgCellType::PRISM:
            if (mapFacePositionToPrism(facePos) == dgFacePositionOnPrism::ABE ||
                mapFacePositionToPrism(facePos) == dgFacePositionOnPrism::CDF)
            {
                JDuffy = ((1 - gaussPt.z())/2)/2;
            }
            break;

        case dgCellType::PYRAMID:
            if (mapFacePositionToPyramid(facePos) == dgFacePositionOnPyramid::ACE ||
                mapFacePositionToPyramid(facePos) == dgFacePositionOnPyramid::CDE ||
                mapFacePositionToPyramid(facePos) == dgFacePositionOnPyramid::BDE ||
                mapFacePositionToPyramid(facePos) == dgFacePositionOnPyramid::ABE)
            {
                JDuffy = ((1 - gaussPt.z())/2)/2;
            }
            break;

        case dgCellType::TET:
            if (mapFacePositionToTet(facePos) == dgFacePositionOnTet::BCE ||
                mapFacePositionToTet(facePos) == dgFacePositionOnTet::ACE ||
                mapFacePositionToTet(facePos) == dgFacePositionOnTet::ABE)
            {
                JDuffy = ((1 - gaussPt.z())/2)/2;
            }
            else if (mapFacePositionToTet(facePos) == dgFacePositionOnTet::ABC)
            {
                JDuffy = ((1 - gaussPt.y())/2)/2;
            }
            break;

        default:
            FatalErrorInFunction
                << "Unsupported cell type for Duffy transformation: "
                << static_cast<int>(cellType) << abort(FatalError);
    }

    return mag(JDuffy);
}
// ************************************************************************* //

