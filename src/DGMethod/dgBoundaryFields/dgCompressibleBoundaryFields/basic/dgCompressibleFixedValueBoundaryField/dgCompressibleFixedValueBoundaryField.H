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

#ifndef Foam_dgCompressibleFixedValueBoundaryField_H
#define Foam_dgCompressibleFixedValueBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class dgCompressibleFixedValueBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleFixedValueBoundaryField
 * \brief Compressible DG boundary condition with a prescribed constant state.
 *
 * \details
 * Reads user-specified primitive data
 * \f$(p_b,\,T_b,\,\mathbf{U}_b)\f$ and converts them once into a
 * conservative ghost state that is reused at every boundary Gauss point.
 *
 * In compact form, the imposed exterior state is
 *
 * \f[
 *     (\rho^+,\,\rho\mathbf{U}^+,\,E^+)
 *     =
 *     \mathcal{C}\!\left(p_b,\,T_b,\,\mathbf{U}_b\right),
 * \f]
 *
 * where \f$\mathcal{C}(\cdot)\f$ denotes the thermo-model conversion from
 * primitive to conservative variables. For an ideal-gas style conservative
 * form this corresponds to
 *
 * \f[
 *     \rho^+ = \rho(p_b, T_b), \qquad
 *     \rho\mathbf{U}^+ = \rho^+ \mathbf{U}_b,
 * \f]
 *
 * together with the compatible total-energy state returned by the selected
 * thermodynamic model.
 */
class dgCompressibleFixedValueBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    /** \brief Prescribed boundary density in conservative form. */
    scalar rhoValue_;

    /** \brief Prescribed boundary momentum in conservative form. */
    vector rhoUValue_;

    /** \brief Prescribed boundary total energy in conservative form. */
    scalar EValue_;

public:

    TypeName("compressible::fixedValue");

    dgCompressibleFixedValueBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressibleFixedValueBoundaryField() = default;

    /**
     * \brief Assign the precomputed constant ghost state.
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
    virtual void updateConservativeGhostState
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

    virtual void updatePrimitiveBCValue
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

    /** \brief Check that the boundary condition is applied on a supported patch type. */
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //