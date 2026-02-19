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
    Foam::dgGeneralInletOutletBoundaryField

Description
    General inlet–outlet boundary condition for DG fields.

SourceFiles
    dgInletOutletBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralInletOutletBoundaryField_H
#define dgGeneralInletOutletBoundaryField_H

#include "dgGeneralBoundaryField.H"
#include "dgField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dgGeneralInletOutletBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeneralInletOutletBoundaryField
 * \brief Inlet–outlet boundary condition for DG fields
 *
 * \details
 * Selects inlet/outlet behaviour per face Gauss point based on the sign of
 * \f$\mathbf{U}\cdot\mathbf{n}\f$ evaluated at that Gauss point.
 *
 * Inlet (\f$\mathbf{U}\cdot\mathbf{n}<0\f$):
 * \f[
 *   \phi^{+} = 2\,\phi_{\mathrm{inlet}} - \phi^{-}
 * \f]
 *
 * Outlet (\f$\mathbf{U}\cdot\mathbf{n}\ge 0\f$):
 * \f[
 *   \phi^{+} = \phi^{-}
 * \f]
 *
 * Requirements:
 * - A velocity field named \c "U" of type dgField<vector> must be available.
 * - The inlet value must be provided via dictionary entry \c inletValue.
 */
template<class Type>
class dgGeneralInletOutletBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
protected:

    /** Prescribed inlet value (Dirichlet) */
    Type inletValue_;

    /** Name of velocity field used for flow-direction check */
    word UName_ = "U";

    /** Velocity field reference (read-only, non-owning) */
    const dgField<vector>& U_;

public:

    //- Runtime type information
    TypeName("generalInletOutlet");

    /**
     * \brief Construct from patch, DG mesh, and dictionary
     *
     * Reads inletValue and resolves the velocity field used for
     * the flow-direction check.
     *
     * \param patch  Finite-volume patch
     * \param dgMesh DG geometric mesh
     * \param dict   Boundary-condition dictionary
     */
    dgGeneralInletOutletBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    /**
     * \brief Update ghost-state value using inlet–outlet logic
     *
     * Uses the sign of \f$\mathbf{U}\cdot\mathbf{n}\f$ at the Gauss point
     * to select between inlet (fixed-value reflection) and outlet
     * (zero-gradient) behaviour.
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
        const vector& Un =
            U_.gaussFields()[cellID].facePlusValueAt(faceLocalID, localGauss);

        const scalar dotProd = Un & n;

        if (dotProd < 0) // Inlet: flow into domain
        {
            plusValue = 2 * inletValue_ - minusValue;
        }
        else // Outlet: flow out of domain
        {
            plusValue = minusValue;
        }
    }

    /**
     * \brief Validate patch compatibility
     *
     * Allowed patch type:
     * - patch
     */
    void checkPatchType() const override
    {
        if (this->patch_.type() != "patch")
        {
            FatalErrorInFunction
                << "Boundary condition " << type() << " can only "
                << "be applied to patch type:\n"
                << "    patch\n"
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

#include "dgGeneralInletOutletBoundaryFieldI.H"

#endif

// ************************************************************************* //
