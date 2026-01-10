/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

Description
    Macros for creating dgInterfaceBoundaryField types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgInterfaceBoundaryFieldMacros_H
#define Foam_dgInterfaceBoundaryFieldMacros_H

#include "dgInterfaceBoundaryField.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#define makedgInterfaceBoundaryField(BCType)                                        \
    defineNamedTemplateTypeNameAndDebug(BCType<scalar>, 0);                      \
    defineNamedTemplateTypeNameAndDebug(BCType<vector>, 0);                      \
    defineNamedTemplateTypeNameAndDebug(BCType<tensor>, 0);                      \
                                                                                 \
    addTemplatedToRunTimeSelectionTable                                          \
    (                                                                            \
        dgInterfaceBoundaryField,                                                \
        BCType,                                                                  \
        scalar,                                                                  \
        fvPatch                                                                  \
    );                                                                           \
                                                                                 \
    addTemplatedToRunTimeSelectionTable                                          \
    (                                                                            \
        dgInterfaceBoundaryField,                                                \
        BCType,                                                                  \
        vector,                                                                  \
        fvPatch                                                                  \
    );                                                                           \
                                                                                 \
    addTemplatedToRunTimeSelectionTable                                          \
    (                                                                            \
        dgInterfaceBoundaryField,                                                \
        BCType,                                                                  \
        tensor,                                                                  \
        fvPatch                                                                  \
    )


#endif

// ************************************************************************* //
