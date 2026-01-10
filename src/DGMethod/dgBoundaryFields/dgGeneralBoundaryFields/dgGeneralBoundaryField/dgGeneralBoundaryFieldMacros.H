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
    Macros for creating dgGeneralBoundaryField types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeneralBoundaryFieldMacros_H
#define Foam_dgGeneralBoundaryFieldMacros_H

#include "dgGeneralBoundaryField.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#define makeDgGeneralBoundaryField(BCType)                                        \
    defineNamedTemplateTypeNameAndDebug(BCType<scalar>, 0);                      \
    defineNamedTemplateTypeNameAndDebug(BCType<vector>, 0);                      \
    defineNamedTemplateTypeNameAndDebug(BCType<tensor>, 0);                      \
                                                                                 \
    addTemplatedToRunTimeSelectionTable                                          \
    (                                                                            \
        dgGeneralBoundaryField,                                                  \
        BCType,                                                                  \
        scalar,                                                                  \
        dictionary                                                               \
    );                                                                           \
                                                                                 \
    addTemplatedToRunTimeSelectionTable                                          \
    (                                                                            \
        dgGeneralBoundaryField,                                                  \
        BCType,                                                                  \
        vector,                                                                  \
        dictionary                                                               \
    );                                                                           \
                                                                                 \
    addTemplatedToRunTimeSelectionTable                                          \
    (                                                                            \
        dgGeneralBoundaryField,                                                  \
        BCType,                                                                  \
        tensor,                                                                  \
        dictionary                                                               \
    )


#endif

// ************************************************************************* //
