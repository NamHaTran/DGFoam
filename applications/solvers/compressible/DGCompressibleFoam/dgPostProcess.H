/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Ha Nam Tran
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

Description
    DG-aware wrapper around OpenFOAM postProcess.H.

    Intended for inclusion inside solver main() before setRootCase.H.
    It reuses the stock functionObject post-processing path, but switches the
    command-line option to -dgPostProcess and constructs the DG mesh/fields
    needed by DGCompressibleFoam-derived function objects.

\*---------------------------------------------------------------------------*/

#ifndef dgPostProcess_H
#define dgPostProcess_H

const Foam::word dgPostProcessPrevOptionName
(
    Foam::argList::postProcessOptionName
);

Foam::argList::postProcessOptionName = "dgPostProcess";

#ifndef CREATE_FIELDS
    #define CREATE_FIELDS createDGMesh.H
    #define DG_POSTPROCESS_UNDEF_CREATE_FIELDS
#endif

#ifndef CREATE_FIELDS_2
    #define CREATE_FIELDS_2 createDGFields.H
    #define DG_POSTPROCESS_UNDEF_CREATE_FIELDS_2
#endif

#ifndef NO_CONTROL
    #define NO_CONTROL
    #define DG_POSTPROCESS_UNDEF_NO_CONTROL
#endif

#include "postProcess.H"

Foam::argList::postProcessOptionName = dgPostProcessPrevOptionName;

#ifdef DG_POSTPROCESS_UNDEF_CREATE_FIELDS
    #undef CREATE_FIELDS
    #undef DG_POSTPROCESS_UNDEF_CREATE_FIELDS
#endif

#ifdef DG_POSTPROCESS_UNDEF_CREATE_FIELDS_2
    #undef CREATE_FIELDS_2
    #undef DG_POSTPROCESS_UNDEF_CREATE_FIELDS_2
#endif

#ifdef DG_POSTPROCESS_UNDEF_NO_CONTROL
    #undef NO_CONTROL
    #undef DG_POSTPROCESS_UNDEF_NO_CONTROL
#endif

#endif

// ************************************************************************* //
