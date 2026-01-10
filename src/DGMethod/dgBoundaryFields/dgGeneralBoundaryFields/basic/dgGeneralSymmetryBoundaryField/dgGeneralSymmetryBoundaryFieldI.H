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

  Class dgGeneralSymmetryBoundaryField Declaration

  Description
    - Implements a general symmetry boundary condition for scalar, vector, and
      tensor fields.

    - This condition assumes the boundary is a symmetry plane:
        * No flux passes through the boundary.
        * Useful for reducing domain size in symmetric problems.

  Usage
    boundaryField
    {
        patchName
        {
            type    generalSymmetry;
        }
    }

  Note
    - Only geometric information (normal vector n) and internal values are
      used. No input parameters required from dictionary.

SourceFiles
    dgGeneralZeroGradientBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralSymmetryBoundaryFieldI_H
#define dgGeneralSymmetryBoundaryFieldI_H

namespace Foam
{

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * //

template<class Type>
dgGeneralSymmetryBoundaryField<Type>::dgGeneralSymmetryBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh,
    const dictionary& dict
)
:
    dgGeneralBoundaryField<Type>(patch, dgMesh, dict)
{}

} // End namespace Foam

#endif

// ************************************************************************* //

