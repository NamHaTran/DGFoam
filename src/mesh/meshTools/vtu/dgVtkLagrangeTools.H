/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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
    Foam::dgVtkLagrange

Description
    Shared native VTK Lagrange export helpers for DG post-processing.

SourceFiles
    dgVtkLagrangeTools.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgVtkLagrangeTools_H
#define Foam_dgVtkLagrangeTools_H

#include "List.H"
#include "dgCellType.H"
#include "fileName.H"
#include "label.H"
#include "point.H"
#include "scalar.H"
#include "vector.H"
#include "word.H"

#include <vector>

namespace Foam
{
namespace dgVtkLagrange
{

struct GridBuffers
{
    std::vector<point> points;
    std::vector<label> connectivity;
    std::vector<label> offsets;
    std::vector<unsigned char> cellTypes;
    std::vector<label> partitions;
};


struct ScalarPointData
{
    word name;
    const std::vector<scalar>* values;
};


struct VectorPointData
{
    word name;
    const std::vector<vector>* values;
};


struct PointDataSchema
{
    word name;
    label nComponents;
};


const char* vtkFloatTypeName();
word cellTypeName(dgCellType type);
bool supportsNativeLagrange(dgCellType type);
unsigned char vtkLagrangeType(dgCellType type);
label exportOrder(label pOrder);
std::vector<label> vtkNodemap(dgCellType type, label nPoints);
std::vector<vector> pyfrStdElementPoints(dgCellType type, label order);
vector xiToEta(dgCellType type, const vector& xi);
void computeBasisAt
(
    dgCellType type,
    const vector& eta,
    label pOrder,
    List<scalar>& basis
);


void writeVtuPiece
(
    const fileName& filePath,
    scalar timeValue,
    const GridBuffers& grid,
    const std::vector<ScalarPointData>& scalarArrays,
    const std::vector<VectorPointData>& vectorArrays
);


void writePvtu
(
    const fileName& filePath,
    scalar timeValue,
    label nProcs,
    const word& pieceStem,
    const std::vector<PointDataSchema>& pointData
);

} // End namespace dgVtkLagrange
} // End namespace Foam

#endif

// ************************************************************************* //
