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
    Foam::dgGeneralBoundaryManager

Description
    Manager for DG general boundary conditions.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeneralBoundaryManager_H
#define Foam_dgGeneralBoundaryManager_H

#include "dgGeomMesh.H"
#include "fvPatch.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "dgGeneralBoundaryField.H"
#include "dgField.H"
#include "dgBoundaryHelper.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgGeneralBoundaryManager Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeneralBoundaryManager
 * \brief Manager for DG general boundary conditions
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * Constructs and stores one dgGeneralBoundaryField per fvPatch and applies
 * them to update ghost (plus-side) values and gradients at face Gauss points.
 *
 * The manager:
 * - Filters patches via dgBoundaryHelper::isGeneralPatch()
 * - Builds boundary-condition objects via runtime selection
 * - Dispatches updateGhostStateValue()/updateGhostStateGrad() per Gauss point
 */
template<class Type>
class dgGeneralBoundaryManager
{
private:

    /** Reference to DG geometric mesh */
    const dgGeomMesh& mesh_;

    /**
     * List of boundary-condition objects (one per fvPatch)
     *
     * The ordering follows fvBoundaryMesh patch ordering.
     */
    DynamicList<autoPtr<dgGeneralBoundaryField<Type>>> bConditions_;

    /** Internal field value used for initialization (uniform) */
    Type internalValue_;

public:

    // Constructors --------------------------------------------------------

    /**
     * \brief Construct from field dictionary and DG mesh
     *
     * \param fieldDict Field dictionary containing boundary definitions
     * \param mesh      DG geometric mesh
     */
    dgGeneralBoundaryManager
    (
        const dictionary& fieldDict,
        const dgGeomMesh& mesh
    );

    /**
     * \brief Construct from IOobject and DG mesh
     *
     * Typically used when reading boundary conditions from disk.
     *
     * \param io   IOobject associated with the field
     * \param mesh DG geometric mesh
     */
    dgGeneralBoundaryManager
    (
        const IOobject& io,
        const dgGeomMesh& mesh
    );

    //- Destructor
    ~dgGeneralBoundaryManager() = default;

    // Access functions ----------------------------------------------------

    /** Return number of boundary patches */
    inline label size() const { return bConditions_.size(); }

    /**
     * \brief Return const access to boundary condition by patch index
     *
     * \param i Patch index
     * \return Const reference to boundary-condition object
     */
    inline const dgGeneralBoundaryField<Type>& operator[](const label i) const
    {
        return bConditions_[i];
    }

    /**
     * \brief Return non-const access to boundary condition by patch index
     *
     * \param i Patch index
     * \return Reference to boundary-condition object
     */
    inline dgGeneralBoundaryField<Type>& operator[](const label i)
    {
        return bConditions_[i];
    }

    /** Return internal field value */
    inline const Type& internalValue() const { return internalValue_; }

    /**
     * \brief Return true if patch is a general (physical) boundary
     *
     * Interface and constraint patches are excluded.
     *
     * \param patchID Patch index
     * \return true if patch is a general boundary
     */
    inline bool isGeneralPatch(const label patchID) const
    {
        const fvBoundaryMesh& patches = mesh_.getFvMesh().boundary();
        return dgBoundaryHelper::isGeneralPatch(patches[patchID]);
    }

    // Core methods --------------------------------------------------------

    /**
     * \brief Update ghost-state values on boundary faces
     *
     * Applies boundary conditions to set plus-side values at face Gauss points.
     *
     * \param GField GaussField to be updated
     */
    void updateValue
    (
        GaussField<Type>& GField
    );

    /**
     * \brief Update ghost-state gradients on boundary faces
     *
     * Applies boundary conditions to set plus-side gradients at face Gauss
     * points.
     *
     * \param GField GaussField to be updated
     */
    void updateGrad
    (
        GaussField<Type>& GField
    );
};

#include "dgGeneralBoundaryManagerI.H"

} // End namespace Foam

#endif
// ************************************************************************* //
