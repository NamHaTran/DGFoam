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

Class
    Foam::dgBasisField

Description
    Stores precomputed GaussField of all basis functions and their derivatives
    for a given cell in DG formulation.

SourceFiles
    dgBasisFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgBasisFieldI_H
#define Foam_dgBasisFieldI_H

namespace Foam
{

// * * * * * * * * Compute scalar basis  * * * * * * * //

inline void dgBasisField::computeBasisField()
{
    const dgGeomCell& cell = *mesh_.cells()[cellID_];
    const List<List<scalar>>& basisCell = cell.basis();

    const label nFaces = faceLabels_.size();
    basis_.setSize(nDof_);

    for (label dof = 0; dof < nDof_; ++dof)
    {
        // Create GaussField for basis[dof]
        GaussField<scalar> phi(cellID_, &mesh_);

        // ---- CELL GAUSS POINTS ---- //
        for (label g = 0; g < cell.nGauss(); ++g)
        {
            phi.cellValueAt(g) = basisCell[g][dof];
        }

        // ---- FACE GAUSS POINTS ---- //
        for (label lf = 0; lf < nFaces; ++lf)
        {
            const label faceID = faceLabels_[lf];
            const dgGeomFace& face = *mesh_.faces()[faceID];
            const bool isBoundary = (faceID >= mesh_.nInternalFaces());
            const bool ownerSide = (mesh_.faceOwner()[faceID] == cellID_);

            // OWNER SIDE CASES
            if (ownerSide)
            {
                const auto& minusB = face.ownerBasis();
                const auto& plusB =
                    isBoundary ? minusB : face.neighborBasis();

                const label nG = face.gaussPointsOwner().size();

                for (label g = 0; g < nG; ++g)
                {
                    phi.faceMinusValueAt(lf, g) = minusB[g][dof];
                    phi.facePlusValueAt(lf, g)  = plusB[g][dof];
                }
            }
            else // NEIGHBOR SIDE CASES
            {
                const auto& minusB  = face.neighborBasis();
                const auto& plusB   = face.ownerBasis();
                const label nG = face.gaussPointsNeighbor().size();

                for (label g = 0; g < nG; ++g)
                {
                    phi.faceMinusValueAt(lf, g) = minusB[g][dof];
                    phi.facePlusValueAt(lf, g)  = plusB[g][dof];
                }
            }
        }

        basis_[dof] = phi;
    }
}


// * * * * * * * Compute derivative basis * * * * * * * //

inline void dgBasisField::computeDBasisField()
{
    const dgGeomCell& cell = *mesh_.cells()[cellID_];

    const auto& dB1 = cell.dBasis_dEta1();
    const auto& dB2 = cell.dBasis_dEta2();
    const auto& dB3 = cell.dBasis_dEta3();

    const label nFaces = faceLabels_.size();
    dBasis_.setSize(nDof_);

    for (label dof = 0; dof < nDof_; ++dof)
    {
        GaussField<vector> dPhi(cellID_, &mesh_);

        // ---- CELL GAUSS POINTS ---- //
        for (label g = 0; g < cell.nGauss(); ++g)
        {
            dPhi.cellValueAt(g) = vector(dB1[g][dof], dB2[g][dof], dB3[g][dof]);
        }

        // ---- FACE GAUSS POINTS ---- //
        for (label lf = 0; lf < nFaces; ++lf)
        {
            const label faceID = faceLabels_[lf];
            const dgGeomFace& face = *mesh_.faces()[faceID];
            const bool isBoundary = (faceID >= mesh_.nInternalFaces());
            const bool ownerSide = (mesh_.faceOwner()[faceID] == cellID_);

            auto loadVec = [&](const List<List<scalar>>& B1,
                               const List<List<scalar>>& B2,
                               const List<List<scalar>>& B3,
                               label g)
            {
                return vector(B1[g][dof], B2[g][dof], B3[g][dof]);
            };

            if (ownerSide)
            {
                const auto& m1 = face.owner_dBasis_dEta1();
                const auto& m2 = face.owner_dBasis_dEta2();
                const auto& m3 = face.owner_dBasis_dEta3();

                const auto& p1 =
                    isBoundary ? m1 : face.neighbor_dBasis_dEta1();
                const auto& p2 =
                    isBoundary ? m2 : face.neighbor_dBasis_dEta2();
                const auto& p3 =
                    isBoundary ? m3 : face.neighbor_dBasis_dEta3();

                const label nG = face.gaussPointsOwner().size();

                for (label g = 0; g < nG; ++g)
                {
                    dPhi.faceMinusValueAt(lf, g) = loadVec(m1, m2, m3, g);
                    dPhi.facePlusValueAt(lf, g)  = loadVec(p1, p2, p3, g);
                }
            }
            else
            {
                const auto& m1 = face.neighbor_dBasis_dEta1();
                const auto& m2 = face.neighbor_dBasis_dEta2();
                const auto& m3 = face.neighbor_dBasis_dEta3();

                const auto& p1 = face.owner_dBasis_dEta1();
                const auto& p2 = face.owner_dBasis_dEta2();
                const auto& p3 = face.owner_dBasis_dEta3();

                const label nG = face.gaussPointsNeighbor().size();

                for (label g = 0; g < nG; ++g)
                {
                    dPhi.faceMinusValueAt(lf, g) = loadVec(m1, m2, m3, g);
                    dPhi.facePlusValueAt(lf, g)  = loadVec(p1, p2, p3, g);
                }
            }
        }

        dBasis_[dof] = dPhi;
    }
}


// Access functions

inline const GaussField<scalar>&
dgBasisField::getBasis(const label dof) const
{
    return basis_[dof];
}

inline const GaussField<vector>&
dgBasisField::getDBasis(const label dof) const
{
    return dBasis_[dof];
}

} // End namespace Foam

#endif

// ************************************************************************* //

