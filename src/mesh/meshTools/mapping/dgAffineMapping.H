/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

Description
    Mapping between reference and physical coordinates for DG affine cells

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgAffineMapping_H
#define Foam_dgAffineMapping_H

#include "fvMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Utility Functions for Affine Mappings
\*---------------------------------------------------------------------------*/

    /**
     * \brief Mapping a point from reference space to physical space for an hexahedral element
     * 
     * \param eta Reference coordinates (eta1, eta2, eta3)
     * \param cellVertices Physical coordinates of the hex vertices (ordered as A, B, C, D, E, F, G, H)
     */
    inline vector mapEtaToXHex(
        const vector& eta, 
        const List<vector>& cellVertices
    );

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
    );

    /**
     * \brief Mapping a point from reference space (eta) to physical space (x) for a prismatic element
     * 
     * \param eta Reference coordinates (eta1, eta2, eta3)
     * \param cellVertices Physical coordinates of the prism vertices (ordered as E, A, B, F, D, C)
     */
    inline vector mapEtaToXPrism(
        const vector& eta, 
        const List<vector>& cellVertices
    );

    /**
     * \brief Map a physical point to reference coordinates of a 6-node prism
     *
     * \details
     * Solves x = sum_i N_i(eta) Xi using Newton iteration.
     * The mapping is trilinear and non-affine.
     *
     * \param x              Physical point
     * \param cellVertices   List of 6 prism vertices [E,A,B,F,D,C]
     *
     * \return Reference coordinates (eta1, eta2, eta3)
     */
    inline vector mapXToEtaPrism
    (
        const vector& x,
        const List<vector>& cellVertices
    );

    /**
     * \brief Mapping a point from reference space (eta) to physical space (x) for a tetrahedral element
     * 
     * \param eta Reference coordinates (eta1, eta2, eta3)
     * \param cellVertices Physical coordinates of the tetrahedron vertices (ordered as A, B, C, D)
     */
    inline vector mapEtaToXTet(
        const vector& eta, 
        const List<vector>& cellVertices
    );

    /**
     * \brief Map a physical point to collapsed reference coordinates of a tetrahedron
     *
     * \details
     * Computes barycentric coordinates by solving the affine mapping:
     *
     *   x = A + λB (B-A) + λC (C-A) + λD (D-A)
     *
     * Then converts barycentric values to collapsed coordinates:
     *
     *   eta3 = 2*λD - 1
     *   eta2 = 2*λC/(1-λD) - 1
     *   eta1 = 2*λB/(1-λC-λD) - 1
     *
     * \param x              Physical point
     * \param cellVertices   List of 4 tetra vertices [A,B,C,D]
     *
     * \return Reference coordinates (eta1, eta2, eta3)
     */
    inline vector mapXToEtaTet
    (
        const vector& x,
        const List<vector>& cellVertices
    );

    /**
     * \brief Mapping a point from reference space (eta) to physical space (x) for a pyramid element
     * 
     * \param eta Reference coordinates (eta1, eta2, eta3)
     * \param cellVertices Physical coordinates of the pyramid vertices (ordered as A, B, C, D, E)
     */
    inline vector mapEtaToXPyramid(
        const vector& eta, 
        const List<vector>& cellVertices
    );

    /**
     * \brief Map a physical point to reference coordinates of a 5-node pyramid
     *
     * \details
     * Solves x = sum_i N_i(eta) Xi using Newton iteration.
     * The mapping is trilinear and non-affine.
     *
     * \param x              Physical point
     * \param cellVertices   List of 5 pyramid vertices [A,B,C,D,E]
     *
     * \return Reference coordinates (eta1, eta2, eta3)
     */
    inline vector mapXToEtaPyramid
    (
        const vector& x,
        const List<vector>& cellVertices
    );

    /**
     * \brief General mapping function that dispatches to the appropriate mapping based on cell type
     * 
     * \param eta Reference coordinates (eta1, eta2, eta3)
     * \param cellVertices Physical coordinates of the cell vertices (ordered according to the cell type)
     * \param cellType Type of the cell (HEXA, PRISM, PYRAMID, TETRA)
     */
    inline vector mapEtaToX(
        const vector& eta, 
        const List<vector>& cellVertices,
        const dgCellType cellType
    );

    /**
     * \brief General mapping function that dispatches to the appropriate mapping based on cell type
     * 
     * \param x Physical coordinates (x1, x2, x3)
     * \param cellVertices Physical coordinates of the cell vertices (ordered according to the cell type)
     * \param cellType Type of the cell (HEXA, PRISM, PYRAMID, TETRA)
     */
    inline vector mapXToEta(
        const vector& x, 
        const List<vector>& cellVertices,
        const dgCellType cellType
    );

    /**
     * \brief Mapping a point from reference space (xi) to reference space (eta, ref hexahedral) for a prismatic element
     */
    inline vector mapXiToEtaPrism(
        const vector& xi 
    );

    /**
     * \brief Mapping a point from reference space (xi) to reference space (eta, ref hexahedral) for a pyramid element
     */
    inline vector mapXiToEtaPyramid(
        const vector& xi 
    );

    /**
     * \brief Mapping a point from reference space (xi) to reference space (eta, ref hexahedral) for a tetrahedral element
     */
    inline vector mapXiToEtaTetra(
        const vector& xi 
    );

    /**
     * \brief Mapping a point from reference space (eta, ref hexahedral) to reference space (xi) for a tetrahedral element
     */
    inline vector mapEtaToXiTetra(
        const vector& eta 
    );

    /**
     * \brief Mapping a point from reference space (eta, ref hexahedral) to reference space (xi) for a prismatic element
     */
    inline vector mapEtaToXiPrism(
        const vector& eta 
    );

    /**
     * \brief Mapping a point from reference space (eta, ref hexahedral) to reference space (xi) for a pyramid element
     */
    inline vector mapEtaToXiPyramid(
        const vector& eta 
    );

    /**
     * \brief Map reference triangle coordinates to physical space
     *
     * \details
     * Computes the physical point X from barycentric coordinates (r, s)
     * on the reference triangle using:
     *
     *   X = (1 - r - s) P0 + r P1 + s P2
     *
     * \param rs  Reference coordinates (r, s) in the reference triangle
     * \param vertices Physical coordinates of the triangle vertices (ordered as A, B, C)
     *
     * \return Physical 3D coordinate
     */
    inline Foam::vector mapToPhysicalTriangle
    (
        const vector2D& rs,
        const List<vector>& vertices
    );

    /**
     * \brief Map reference square coordinates to physical quadrilateral
     *
     * \details
     * Uses bilinear shape functions for a 4-node quad.
     *
     * \param xiEta  Reference coordinates (xi, eta) in [-1,1]^2
     * \param vertices Physical coordinates of the quadrilateral vertices (ordered as A, B, C, D)
     *
     * \return Physical 3D coordinate
     */
    inline Foam::vector mapToPhysicalQuad
    (
        const vector2D& eta,
        const List<vector>& vertices
    );

} // End namespace Foam

#include "dgAffineMappingI.H"

#endif

// ************************************************************************* //
