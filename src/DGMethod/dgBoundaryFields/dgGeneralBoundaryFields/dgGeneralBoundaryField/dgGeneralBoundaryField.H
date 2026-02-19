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
    Foam::dgGeneralBoundaryField

Description
    Base class for runtime-selectable DG boundary conditions.

SourceFiles
    dgGeneralBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeneralBoundaryField_H
#define Foam_dgGeneralBoundaryField_H

#include "word.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dgGeomMesh.H"
#include "fvPatch.H"

namespace Foam
{

// Forward declaration
template<class Type> class faceGaussField;

/*---------------------------------------------------------------------------*\
                    Class dgGeneralBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeneralBoundaryField
 * \brief Base class for DG ghost-state boundary conditions
 *
 * \details
 * Defines a runtime-selectable interface for DG boundary conditions that
 * compute the exterior (plus-side) state from the interior (minus-side)
 * state at boundary face Gauss points:
 *
 * \f[
 *   \phi^{-} \;\longrightarrow\; \phi^{+}
 * \f]
 *
 * Derived classes typically override:
 * - updateGhostStateValue() to define the ghost-state value
 * - updateGhostStateGrad()  to define the ghost-state gradient
 *
 * The default implementation corresponds to a copy/zero-gradient behaviour.
 */
template<class Type>
class dgGeneralBoundaryField
{
protected:

    /** Reference to fvPatch describing boundary geometry */
    const fvPatch& patch_;

    /** Reference to DG geometric mesh */
    const dgGeomMesh& dgMesh_;

    /** Boundary-condition dictionary (stored by value) */
    const dictionary dict_;

public:

    //- Runtime type information
    TypeName("dgGeneralBoundaryField");

    /**
     * \brief Runtime selection table for DG boundary conditions
     *
     * Enables dictionary-based construction of derived boundary-condition
     * types.
     */
    declareRunTimeSelectionTable
    (
        autoPtr,
        dgGeneralBoundaryField,
        dictionary,
        (const fvPatch& patch, const dgGeomMesh& dgMesh, const dictionary& dict),
        (patch, dgMesh, dict)
    );

    // Constructors --------------------------------------------------------

    /**
     * \brief Construct from patch, DG mesh, and dictionary
     *
     * \param patch  Finite-volume patch
     * \param dgMesh DG geometric mesh
     * \param dict   Boundary-condition dictionary
     */
    dgGeneralBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralBoundaryField() = default;

    // Factory method ------------------------------------------------------

    /**
     * \brief Create boundary condition from dictionary
     *
     * Uses the runtime-selection mechanism to construct a derived
     * boundary-condition type.
     *
     * \param patch  Finite-volume patch
     * \param dgMesh DG geometric mesh
     * \param dict   Boundary-condition dictionary
     * \return autoPtr to constructed dgGeneralBoundaryField<Type>
     */
    static autoPtr<dgGeneralBoundaryField<Type>> New
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    // Core DG boundary interface -----------------------------------------

    /**
     * \brief Update ghost-state value for scalar fields
     *
     * Default behaviour copies the interior value.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusValue    Interior (minus-side) value
     * \param plusValue     Ghost (plus-side) value to be set
     */
    virtual void updateGhostStateValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar& minusValue,
        scalar& plusValue
    ) const
    {
        plusValue = minusValue;
    }

    /**
     * \brief Update ghost-state value for vector fields
     *
     * Default behaviour copies the interior value.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusValue    Interior (minus-side) value
     * \param plusValue     Ghost (plus-side) value to be set
     */
    virtual void updateGhostStateValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const vector& minusValue,
        vector& plusValue
    ) const
    {
        plusValue = minusValue;
    }

    /**
     * \brief Update ghost-state value for tensor fields
     *
     * Default behaviour copies the interior value.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusValue    Interior (minus-side) value
     * \param plusValue     Ghost (plus-side) value to be set
     */
    virtual void updateGhostStateValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const tensor& minusValue,
        tensor& plusValue
    ) const
    {
        plusValue = minusValue;
    }

    /**
     * \brief Update ghost-state gradient for scalar fields
     *
     * Default behaviour copies the interior gradient.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusGrad     Interior (minus-side) gradient
     * \param plusGrad      Ghost (plus-side) gradient to be set
     */
    virtual void updateGhostStateGrad
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar& minusGrad,
        scalar& plusGrad
    ) const
    {
        plusGrad = minusGrad;
    }

    /**
     * \brief Update ghost-state gradient for vector fields
     *
     * Default behaviour copies the interior gradient.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusGrad     Interior (minus-side) gradient
     * \param plusGrad      Ghost (plus-side) gradient to be set
     */
    virtual void updateGhostStateGrad
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const vector& minusGrad,
        vector& plusGrad
    ) const
    {
        plusGrad = minusGrad;
    }

    /**
     * \brief Update ghost-state gradient for tensor fields
     *
     * Default behaviour copies the interior gradient.
     *
     * \param cellID        Global cell index
     * \param faceLocalID   Local face index within the cell
     * \param localGauss    Local Gauss-point index on the face
     * \param n             Outward unit normal
     * \param minusGrad     Interior (minus-side) gradient
     * \param plusGrad      Ghost (plus-side) gradient to be set
     */
    virtual void updateGhostStateGrad
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const tensor& minusGrad,
        tensor& plusGrad
    ) const
    {
        plusGrad = minusGrad;
    }

    // Access functions ----------------------------------------------------

    /** Return boundary patch */
    inline const fvPatch& patch() const { return patch_; }

    /** Return DG geometric mesh */
    inline const dgGeomMesh& dgMesh() const { return dgMesh_; }

    /** Return boundary-condition dictionary */
    inline const dictionary& dict() const { return dict_; }

    /** Return patch name */
    inline const word& name() const { return patch_.name(); }

    /**
     * \brief Return boundary-condition type name
     *
     * Extracts the value of dictionary entry \c type.
     *
     * \return Boundary-condition type name
     */
    inline word BCtype() const
    {
        word t;
        dict_.lookup("type") >> t;
        return t;
    }

    // Patch compatibility -------------------------------------------------

    /**
     * \brief Validate patch compatibility
     *
     * Derived classes should override this method to restrict supported
     * patch types.
     */
    virtual void checkPatchType() const {}
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgGeneralBoundaryFieldI.H"

#endif

// ************************************************************************* //
