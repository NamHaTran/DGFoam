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

Description
    Tools for dgGeomFace

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeomFaceTools_H
#define Foam_dgGeomFaceTools_H

#include "fvMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Utility Functions for Face–Patch Mapping
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
    Description
        Return the patch ID that owns the given faceID.

        - If faceID belongs to an internal face (faceID < nInternalFaces),
          return -1.
        - If boundary face, detect which fvPatch owns it using the fact that
          fvPatch stores a consecutive block of boundary faces:
              [start, start + size)

        This is the same logic OpenFOAM uses internally for patch ownership.

    Arguments
        mesh   : reference fvMesh providing boundary and internal face layout
        faceID : global face ID in fvMesh.faces()

    Return
        label : patch index, or -1 for internal face
\*---------------------------------------------------------------------------*/
inline label findOwnerPatch(const fvMesh& mesh, const label faceID);

} // End namespace Foam

#include "dgGeomFaceToolsI.H"

#endif

// ************************************************************************* //
