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

Description
    Helper utilities for classification of boundary patches in DG.

SourceFiles
    dgBoundaryHelper.H

\*---------------------------------------------------------------------------*/

#ifndef dgBoundaryHelper_H
#define dgBoundaryHelper_H

#include "fvPatch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace dgBoundaryHelper
\*---------------------------------------------------------------------------*/

/**
 * \namespace Foam::dgBoundaryHelper
 * \brief Utility functions for DG boundary patch classification.
 *
 * \details
 * This namespace provides a lightweight set of helper functions
 * to classify OpenFOAM boundary patches into logical categories
 * used by the DG framework:
 *  - Interface patches (e.g. processor boundaries)
 *  - General physical boundary patches (walls, inlets, outlets, symmetry)
 *  - Constraint or special patches (e.g. overset, empty)
 *
 * The classification is based solely on the fvPatch type string and is
 * intentionally kept simple and extensible.
 */
namespace dgBoundaryHelper
{

//======================================================================//
//  Interface boundary patches
//======================================================================//

/**
 * \brief Return true if the patch is treated as an interface patch.
 *
 * Interface patches are used to connect DG cells across
 * processor boundaries or other internal interfaces.
 *
 * \param patch Boundary patch to classify.
 * \return True for DG interface patches.
 */
inline bool isInterfacePatch(const fvPatch& patch)
{
    const word& type = patch.type();

    return
    (
        type == "processor"
    );
}


//======================================================================//
//  General (physical) boundary patches
//======================================================================//

/**
 * \brief Return true if the patch is treated as a general physical boundary.
 *
 * General patches represent physical boundaries such as
 * walls, inlets, outlets, and symmetry planes.
 *
 * \param patch Boundary patch to classify.
 * \return True for general physical boundary patches.
 */
inline bool isGeneralPatch(const fvPatch& patch)
{
    const word& type = patch.type();

    return
    (
        type == "wall"
     || type == "patch"          // generic patch (inlet, outlet, etc.)
     || type == "symmetry"
     || type == "symmetryPlane"
    );
}


//======================================================================//
//  Constraint / special boundary patches
//======================================================================//

/**
 * \brief Return true if the patch is treated as a constraint or special patch.
 *
 * Constraint patches impose additional topological or
 * numerical constraints, or mark directions that should not receive
 * DG ghost-state updates (for example, overset and empty).
 *
 * \param patch Boundary patch to classify.
 * \return True for constraint or special patches.
 */
inline bool isConstraintPatch(const fvPatch& patch)
{
    const word& type = patch.type();

    return
    (
        type == "overset"
     || type == "empty"
    );
}


//======================================================================//
//  Sanity check helpers
//======================================================================//

/**
 * \brief Return true if the patch type is recognised by the DG framework.
 *
 * Returns true if the patch belongs to any known DG category.
 *
 * \param patch Boundary patch to classify.
 * \return True if the patch belongs to a supported DG category.
 */
inline bool isKnownPatchType(const fvPatch& patch)
{
    return
    (
        isGeneralPatch(patch)
     || isInterfacePatch(patch)
     || isConstraintPatch(patch)
    );
}

} // End namespace dgBoundaryHelper

} // End namespace Foam

#endif

// ************************************************************************* //
