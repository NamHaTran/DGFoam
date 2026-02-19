/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeomFaceToolsI_H
#define Foam_dgGeomFaceToolsI_H

namespace Foam
{

inline label findOwnerPatch(const fvMesh& mesh, const label faceID)
{
    const label nInternal = mesh.nInternalFaces();

    // Internal faces never belong to a patch
    if (faceID < nInternal)
    {
        return -1;
    }

    const fvBoundaryMesh& patches = mesh.boundary();

    // Loop over all boundary patches and test patch range
    for (label patchI = 0; patchI < patches.size(); ++patchI)
    {
        const fvPatch& p = patches[patchI];
        const label start = p.start();
        const label end   = start + p.size();

        if (faceID >= start && faceID < end)
        {
            return patchI;
        }
    }

    // Should not happen unless mesh is corrupted
    return -1;
}

} // End namespace Foam

#endif

// ************************************************************************* //
