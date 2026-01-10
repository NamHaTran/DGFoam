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
    Helper functions for boundary patches classification.

\*---------------------------------------------------------------------------*/

#ifndef dgBoundaryHelper_H
#define dgBoundaryHelper_H

#include "fvPatch.H"

namespace Foam
{

namespace dgBoundaryHelper
{

//======================================================================//
//  Interface boundary patches
//======================================================================//

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

inline bool isGeneralPatch(const fvPatch& patch)
{
    const word& type = patch.type();

    return
    (
        type == "wall"
     || type == "patch"        // generic patch (inlet, outlet, etc.)
     || type == "symmetry"
     || type == "symmetryPlane"
    );
}


//======================================================================//
//  Constraint / special boundary patches (future extension)
//======================================================================//

inline bool isConstraintPatch(const fvPatch& patch)
{
    const word& type = patch.type();

    return
    (
        type == "overset"
    );
}


//======================================================================//
//  Sanity check helpers
//======================================================================//

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
