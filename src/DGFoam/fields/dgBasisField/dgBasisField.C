/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "dgBasisField.H"
#include "dgGeomMesh.H"
#include "dgGeomCell.H"
#include "dgGeomFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::dgBasisField::dgBasisField
(
    const label cellID,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    cellID_(cellID),
    nDof_(mesh_.cells()[cellID_]->nDof())
{
    // Store local face labels for this cell
    faceLabels_ = mesh_.cells()[cellID_]->faces();

    // Precompute basis functions at Gauss points
    computeBasisField();
    computeDBasisField();
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::dgBasisField::computeBasisField()
{
    const dgGeomCell& cell = *mesh_.cells()[cellID_];
    const List<List<scalar>>& basis = cell.basis();

    //const label nDof = cell.nDof();
    const label nFaces = faceLabels_.size();

    basis_.setSize(nDof_);

    for (label dof = 0; dof < nDof_; ++dof)
    {
        // Initialize GaussField for current basis function
        GaussField<scalar> phi(cellID_, &mesh_);

        // --- Cell Gauss points --- //
        for (label nG = 0; nG < cell.nGauss(); ++nG)
        {
            phi.cellValueAt(nG) = basis[nG][dof];
        }

        // --- Face Gauss points --- //
        for (label faceI = 0; faceI < nFaces; ++faceI)
        {
            const label faceID = faceLabels_[faceI];
            const dgGeomFace& face = *mesh_.faces()[faceID];
            const label owner = mesh_.faceOwner()[faceID];
            const bool isBoundaryFace = (faceID >= mesh_.nInternalFaces());

            if (owner == cellID_)
            {
                if (isBoundaryFace)
                {
                    const List<List<scalar>>& basisMinus = face.ownerBasis();
                    const List<List<scalar>>& basisPlus = basisMinus;
                    const label nGauss = face.gaussPointsOwner().size();

                    for (label nG = 0; nG < nGauss; ++nG)
                    {
                        phi.faceMinusValueAt(faceI, nG) = basisMinus[nG][dof];
                        phi.facePlusValueAt(faceI, nG)  = basisPlus[nG][dof];
                    }
                }
                else
                {
                    const List<List<scalar>>& basisMinus = face.ownerBasis();
                    const List<List<scalar>>& basisPlus = face.neighborBasis();
                    const label nGauss = face.gaussPointsOwner().size();

                    for (label nG = 0; nG < nGauss; ++nG)
                    {
                        phi.faceMinusValueAt(faceI, nG) = basisMinus[nG][dof];
                        phi.facePlusValueAt(faceI, nG)  = basisPlus[nG][dof];
                    }
                }
            }
            else
            {
                const List<List<scalar>>& basisMinus = face.neighborBasis();
                const List<List<scalar>>& basisPlus  = face.ownerBasis();

                const label nGauss = face.gaussPointsNeighbor().size();

                for (label nG = 0; nG < nGauss; ++nG)
                {
                    phi.faceMinusValueAt(faceI, nG) = basisMinus[nG][dof];
                    phi.facePlusValueAt(faceI, nG)  = basisPlus[nG][dof];
                }
            }
        }

        basis_[dof] = phi;
    }
}


void Foam::dgBasisField::computeDBasisField()
{
    const dgGeomCell& cell = *mesh_.cells()[cellID_];

    const List<List<scalar>>& dB1 = cell.dBasis_dEta1();
    const List<List<scalar>>& dB2 = cell.dBasis_dEta2();
    const List<List<scalar>>& dB3 = cell.dBasis_dEta3();

    //const label nDof = cell.nDof();
    const label nFaces = faceLabels_.size();

    dBasis_.setSize(nDof_);

    for (label dof = 0; dof < nDof_; ++dof)
    {
        // Initialize GaussField for current derivative basis
        GaussField<vector> dPhi(cellID_, &mesh_);

        // --- Cell Gauss points --- //
        for (label nG = 0; nG < cell.nGauss(); ++nG)
        {
            dPhi.cellValueAt(nG) = vector
            (
                dB1[nG][dof],
                dB2[nG][dof],
                dB3[nG][dof]
            );
        }

        // --- Face Gauss points --- //
        for (label faceI = 0; faceI < nFaces; ++faceI)
        {
            const label faceID = faceLabels_[faceI];
            const dgGeomFace& face = *mesh_.faces()[faceID];
            const label owner = mesh_.faceOwner()[faceID];
            const bool isBoundaryFace = (faceID >= mesh_.nInternalFaces());

            if (owner == cellID_)
            {
                

                if (isBoundaryFace)
                {
                    const List<List<scalar>>& dB1minus = face.owner_dBasis_dEta1();
                    const List<List<scalar>>& dB2minus = face.owner_dBasis_dEta2();
                    const List<List<scalar>>& dB3minus = face.owner_dBasis_dEta3();

                    const List<List<scalar>>& dB1plus  = dB1minus;
                    const List<List<scalar>>& dB2plus  = dB2minus;
                    const List<List<scalar>>& dB3plus  = dB3minus;

                    const label nGauss = face.gaussPointsOwner().size();

                    for (label nG = 0; nG < nGauss; ++nG)
                    {
                        const vector dMinus
                        (
                            dB1minus[nG][dof],
                            dB2minus[nG][dof],
                            dB3minus[nG][dof]
                        );

                        const vector dPlus
                        (
                            dB1plus[nG][dof],
                            dB2plus[nG][dof],
                            dB3plus[nG][dof]
                        );

                        dPhi.faceMinusValueAt(faceI, nG) = dMinus;
                        dPhi.facePlusValueAt(faceI, nG)  = dPlus;
                    }
                }
                else
                {
                    const List<List<scalar>>& dB1minus = face.owner_dBasis_dEta1();
                    const List<List<scalar>>& dB2minus = face.owner_dBasis_dEta2();
                    const List<List<scalar>>& dB3minus = face.owner_dBasis_dEta3();

                    const List<List<scalar>>& dB1plus  = face.neighbor_dBasis_dEta1();
                    const List<List<scalar>>& dB2plus  = face.neighbor_dBasis_dEta2();
                    const List<List<scalar>>& dB3plus  = face.neighbor_dBasis_dEta3();

                    const label nGauss = face.gaussPointsOwner().size();

                    for (label nG = 0; nG < nGauss; ++nG)
                    {
                        const vector dMinus
                        (
                            dB1minus[nG][dof],
                            dB2minus[nG][dof],
                            dB3minus[nG][dof]
                        );

                        const vector dPlus
                        (
                            dB1plus[nG][dof],
                            dB2plus[nG][dof],
                            dB3plus[nG][dof]
                        );

                        dPhi.faceMinusValueAt(faceI, nG) = dMinus;
                        dPhi.facePlusValueAt(faceI, nG)  = dPlus;
                    }
                }
            }
            else
            {
                const List<List<scalar>>& dB1minus = face.neighbor_dBasis_dEta1();
                const List<List<scalar>>& dB2minus = face.neighbor_dBasis_dEta2();
                const List<List<scalar>>& dB3minus = face.neighbor_dBasis_dEta3();

                const List<List<scalar>>& dB1plus  = face.owner_dBasis_dEta1();
                const List<List<scalar>>& dB2plus  = face.owner_dBasis_dEta2();
                const List<List<scalar>>& dB3plus  = face.owner_dBasis_dEta3();

                const label nGauss = face.gaussPointsNeighbor().size();

                for (label nG = 0; nG < nGauss; ++nG)
                {
                    const vector dMinus
                    (
                        dB1minus[nG][dof],
                        dB2minus[nG][dof],
                        dB3minus[nG][dof]
                    );

                    const vector dPlus
                    (
                        dB1plus[nG][dof],
                        dB2plus[nG][dof],
                        dB3plus[nG][dof]
                    );

                    dPhi.faceMinusValueAt(faceI, nG) = dMinus;
                    dPhi.facePlusValueAt(faceI, nG)  = dPlus;
                }
            }
        }

        dBasis_[dof] = dPhi;
    }
}


const Foam::GaussField<scalar>& Foam::dgBasisField::getBasis(const label dof) const
{
    return basis_[dof];
}


const Foam::GaussField<vector>& Foam::dgBasisField::getDBasis(const label dof) const
{
    return dBasis_[dof];
}


// ************************************************************************* //

