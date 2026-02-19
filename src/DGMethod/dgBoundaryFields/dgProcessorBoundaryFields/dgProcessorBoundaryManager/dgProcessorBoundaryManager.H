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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgProcessorBoundaryManager

Description
    Manager for DG processor boundary communication.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgProcessorBoundaryManager_H
#define Foam_dgProcessorBoundaryManager_H

#include "dgGeomMesh.H"
#include "fvPatch.H"
#include "PtrList.H"
#include "dgProcessorBoundaryField.H"
#include "dgField.H"
#include "dgBoundaryHelper.H"
#include "Pstream.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dgProcessorBoundaryManager Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgProcessorBoundaryManager
 * \brief Manager for DG processor boundary fields
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * Detects processor patches in the fvMesh and constructs one
 * dgProcessorBoundaryField per processor patch.
 *
 * The manager coordinates MPI exchanges to synchronise plus-side
 * Gauss values on processor faces. Physical boundary conditions are
 * handled separately by general boundary managers.
 */
template<class Type>
class dgProcessorBoundaryManager
{
private:

    /** Reference to DG geometric mesh */
    const dgGeomMesh& mesh_;

    /**
     * List of processor boundary handlers
     *
     * One dgProcessorBoundaryField is created per processor fvPatch.
     * Ordering follows fvBoundaryMesh patch ordering.
     */
    PtrList<dgProcessorBoundaryField<Type>> bConditions_;

public:

    // Constructors --------------------------------------------------------

    /**
     * \brief Construct processor boundary manager from DG mesh
     *
     * Detects all processor patches from the underlying fvMesh.
     *
     * \param mesh DG geometric mesh
     */
    dgProcessorBoundaryManager
    (
        const dgGeomMesh& mesh
    );

    //- Destructor
    ~dgProcessorBoundaryManager() = default;

    // Access --------------------------------------------------------------

    /** Return number of processor boundary patches */
    inline label size() const
    {
        return bConditions_.size();
    }

    /**
     * \brief Return const access to processor boundary field by index
     *
     * \param i Patch index
     * \return Const reference to processor boundary handler
     */
    inline const dgProcessorBoundaryField<Type>& operator[](const label i) const
    {
        return bConditions_[i];
    }

    /**
     * \brief Return non-const access to processor boundary field by index
     *
     * \param i Patch index
     * \return Reference to processor boundary handler
     */
    inline dgProcessorBoundaryField<Type>& operator[](const label i)
    {
        return bConditions_[i];
    }

    // Core functionality -------------------------------------------------

    /**
     * \brief Update DG field across all processor boundaries
     *
     * Triggers communication on all processor patches to synchronise
     * plus-side Gauss values.
     *
     * \param field DG field to be updated
     */
    void update
    (
        dgField<Type>& field
    );
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgProcessorBoundaryManagerI.H"

#endif
// ************************************************************************* //
