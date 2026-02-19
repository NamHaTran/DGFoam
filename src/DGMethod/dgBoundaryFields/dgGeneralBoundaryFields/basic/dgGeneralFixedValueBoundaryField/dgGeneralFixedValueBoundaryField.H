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
    Foam::dgGeneralFixedValueBoundaryField

Description
    Fixed-value (Dirichlet) boundary condition for DG methods.

SourceFiles
    dgGeneralFixedValueBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralFixedValueBoundaryField_H
#define dgGeneralFixedValueBoundaryField_H

#include "dgGeneralBoundaryField.H"
#include "wallPolyPatch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dgGeneralFixedValueBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeneralFixedValueBoundaryField
 * \brief Strong Dirichlet boundary condition for DG fields
 *
 * \details
 * Enforces a prescribed boundary value using a ghost-state reflection.
 *
 * The Dirichlet condition is imposed strongly by mirroring the
 * interior (minus-side) state about the prescribed value:
 *
 * \f[
 *   U^{+} = 2\,U_{\mathrm{bc}} - U^{-}
 * \f]
 *
 * where:
 * - \f$U^{-}\f$ is the interior (minus-side) value
 * - \f$U^{+}\f$ is the ghost (plus-side) value
 * - \f$U_{\mathrm{bc}}\f$ is the prescribed boundary value
 *
 * No explicit gradient modification is applied; the boundary value
 * fully constrains the solution at the interface.
 */
template<class Type>
class dgGeneralFixedValueBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
    /** Prescribed boundary value */
    Type value_;

public:

    //- Runtime type information
    TypeName("generalFixedValue");

    /**
     * \brief Construct from patch, DG mesh, and dictionary
     *
     * Reads the prescribed boundary value from the dictionary and
     * validates the patch type.
     *
     * \param patch  Finite-volume patch
     * \param dgMesh DG geometric mesh
     * \param dict   Boundary-condition dictionary
     */
    dgGeneralFixedValueBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralFixedValueBoundaryField() = default;

    /**
     * \brief Update ghost-state value for fixed-value condition
     *
     * Applies the strong Dirichlet relation:
     * \f$ U^{+} = 2\,U_{\mathrm{bc}} - U^{-} \f$
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward face normal
     * \param minusValue    Interior (minus-side) value
     * \param plusValue     Ghost (plus-side) value to be set
     */
    virtual void updateGhostStateValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const Type& minusValue,
        Type& plusValue
    ) const override
    {
        plusValue = 2 * value_ - minusValue;
    }

    /**
     * \brief Validate patch compatibility
     *
     * Allowed patch types:
     * - patch
     * - wall
     */
    void checkPatchType() const override
    {
        if
        (
            this->patch_.type() != "patch"
         && this->patch_.type() != "wall"
        )
        {
            FatalErrorInFunction
                << "Boundary condition " << type() << " can only "
                << "be applied to patch types:\n"
                << "    patch\n"
                << "    wall\n"
                << "but patch "
                << this->patch_.name()
                << " is of type " << this->patch_.type() << nl
                << exit(FatalError);
        }
    }
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
              Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgGeneralFixedValueBoundaryFieldI.H"

#endif

// ************************************************************************* //
