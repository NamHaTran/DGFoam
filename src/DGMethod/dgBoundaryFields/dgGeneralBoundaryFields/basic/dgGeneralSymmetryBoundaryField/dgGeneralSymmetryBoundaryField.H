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
    Foam::dgGeneralSymmetryBoundaryField

Description
    Symmetry boundary condition for DG fields.

SourceFiles
    dgGeneralSymmetryBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralSymmetryBoundaryField_H
#define dgGeneralSymmetryBoundaryField_H

#include "dgGeneralBoundaryField.H"
#include "symmetryPolyPatch.H"
#include "vector.H"
#include "tensor.H"
#include "VectorSpace.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dgGeneralSymmetryBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeneralSymmetryBoundaryField
 * \brief Symmetry boundary condition for DG fields
 *
 * \details
 * Enforces a symmetry-plane constraint at face Gauss points by
 * reflecting the interior (minus-side) state with respect to the
 * boundary unit normal.
 *
 * Scalar fields:
 * \f[
 *   \phi^{+} = \phi^{-}
 * \f]
 *
 * Vector fields:
 * \f[
 *   \boldsymbol{\phi}^{+}
 *   =
 *   \boldsymbol{\phi}^{-}
 *   - 2(\boldsymbol{\phi}^{-}\cdot\mathbf{n})\,\mathbf{n}
 * \f]
 *
 * Only interior values and geometric information are used.
 * No dictionary parameters are required.
 */
template<class Type>
class dgGeneralSymmetryBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
public:

    //- Runtime type information
    TypeName("generalSymmetry");

    /**
     * \brief Construct symmetry boundary field
     *
     * \param patch  Finite-volume patch
     * \param dgMesh DG geometric mesh
     * \param dict   Boundary-condition dictionary (unused)
     */
    dgGeneralSymmetryBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralSymmetryBoundaryField() = default;

    /**
     * \brief Update ghost-state value for scalar fields
     *
     * Copies the interior value to the ghost state.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusValue    Interior (minus-side) value
     * \param plusValue     Ghost (plus-side) value to be set
     */
    inline void updateGhostStateValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar& minusValue,
        scalar& plusValue
    ) const override
    {
        plusValue = minusValue;
    }

    /**
     * \brief Update ghost-state value for vector fields
     *
     * Reflects the normal component and preserves tangential components.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusValue    Interior (minus-side) vector value
     * \param plusValue     Ghost (plus-side) vector value to be set
     */
    void updateGhostStateValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const vector& minusValue,
        vector& plusValue
    ) const override
    {
        plusValue = minusValue - 2 * (n & minusValue) * n;
    }

    /**
     * \brief Validate patch compatibility
     *
     * Allowed patch type:
     * - symmetry
     */
    void checkPatchType() const override
    {
        if (this->patch_.type() != "symmetry")
        {
            FatalErrorInFunction
                << "Boundary condition " << type() << " can only "
                << "be applied to patch type:\n"
                << "    symmetry\n"
                << "but patch "
                << this->patch_.name() << " is of type "
                << this->patch_.type() << nl
                << exit(FatalError);
        }
    }
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgGeneralSymmetryBoundaryFieldI.H"

#endif

// ************************************************************************* //
