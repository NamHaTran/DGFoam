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

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCompressibleSymmetryBoundaryField_H
#define Foam_dgCompressibleSymmetryBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class dgCompressibleSymmetryBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleSymmetryBoundaryField
 * \brief Symmetry boundary condition for compressible DG states.
 *
 * \details
 * This boundary condition reflects only the normal component of momentum and
 * leaves density and total energy unchanged:
 *
 * \f[
 *     \rho^+ = \rho^-,
 *     \qquad
 *     E^+ = E^-,
 * \f]
 *
 * \f[
 *     \rho\mathbf{U}^+
 *     =
 *     \rho\mathbf{U}^-
 *     -
 *     2\left(\mathbf{n}\cdot\rho\mathbf{U}^-\right)\mathbf{n}.
 * \f]
 *
 * Equivalently, in primitive variables the tangential velocity is preserved
 * and the normal velocity changes sign across the symmetry plane.
 */
class dgCompressibleSymmetryBoundaryField
:
    public dgCompressibleBoundaryField
{
public:

    TypeName("compressible::symmetry");

    dgCompressibleSymmetryBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressibleSymmetryBoundaryField() = default;

    /**
     * \brief Reflect the interior state across the boundary normal.
     *
     * \param cellID Owner cell index.
     * \param faceLocalID Local face index on the owner cell.
     * \param localGauss Local Gauss-point index on the face.
     * \param n Outward unit normal of the owner cell.
     * \param rhoMinus Interior density.
     * \param rhoUMinus Interior momentum.
     * \param EMinus Interior total energy.
     * \param rhoPlus Ghost density.
     * \param rhoUPlus Ghost momentum.
     * \param EPlus Ghost total energy.
     */
    virtual void updateGhostState
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        scalar& rhoPlus,
        vector& rhoUPlus,
        scalar& EPlus
    ) const override;

    virtual void updateBCValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        scalar& rhoBC,
        vector& rhoUBC,
        scalar& EBC
    ) const override;

    /** \brief Check that the boundary condition is applied on a symmetry patch. */
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
