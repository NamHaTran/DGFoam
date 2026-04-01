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
    Foam::dgEquiSpacedTessellation

Description
    Cached equi-spaced sampling and sub-cell tessellation for DG export.

SourceFiles
    dgEquiSpacedTessellation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgEquiSpacedTessellation_H
#define Foam_dgEquiSpacedTessellation_H

#include "DynamicList.H"
#include "FixedList.H"
#include "List.H"
#include "cellShape.H"
#include "dgCellType.H"
#include "label.H"
#include "pointField.H"
#include "scalar.H"
#include "vector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgEquiSpacedTessellation Declaration
\*---------------------------------------------------------------------------*/

class dgEquiSpacedTessellation
{
public:

    struct SamplePoint
    {
        vector eta;
        List<scalar> basis;
    };

    struct CellInfo
    {
        dgCellType type;
        pointField vertices;
    };

    struct VtkSubCell
    {
        FixedList<label, 8> pointLabels;
        label nPoints;
        unsigned char vtkType;
    };

private:

    label n_;
    label pOrder_;

    mutable List<SamplePoint> hexSamples_;
    mutable List<SamplePoint> prismSamples_;
    mutable List<SamplePoint> pyramidSamples_;
    mutable List<SamplePoint> tetSamples_;
    mutable FixedList<List<SamplePoint>, 3> hexCollapsedSamples_;
    mutable List<SamplePoint> prismCollapsedSamples_;

    mutable List<VtkSubCell> hexSubCells_;
    mutable List<VtkSubCell> prismSubCells_;
    mutable List<VtkSubCell> pyramidSubCells_;
    mutable List<VtkSubCell> tetSubCells_;
    mutable FixedList<List<VtkSubCell>, 3> hexCollapsedSubCells_;
    mutable List<VtkSubCell> prismCollapsedSubCells_;

    static void computeBasisAt
    (
        const dgCellType type,
        const vector& eta,
        const label pOrder,
        List<scalar>& basis
    );

    static vector canonicalEtaFromXi(const vector& xi, const dgCellType type);

    SamplePoint makeSample(const vector& eta, const dgCellType type) const;

    const List<VtkSubCell>& subCells(const dgCellType type) const;
    const List<VtkSubCell>& collapsedSubCells
    (
        const dgCellType type,
        const label collapseDir
    ) const;

    void buildHexSamples(List<SamplePoint>& samples) const;
    void buildPrismSamples(List<SamplePoint>& samples) const;
    void buildPyramidSamples(List<SamplePoint>& samples) const;
    void buildTetSamples(List<SamplePoint>& samples) const;
    void buildHexCollapsedSamples
    (
        List<SamplePoint>& samples,
        const label collapseDir
    ) const;
    void buildPrismCollapsedSamples(List<SamplePoint>& samples) const;

    void buildHexSubCells(List<VtkSubCell>& cells) const;
    void buildPrismSubCells(List<VtkSubCell>& cells) const;
    void buildPyramidSubCells(List<VtkSubCell>& cells) const;
    void buildTetSubCells(List<VtkSubCell>& cells) const;
    void buildHexCollapsedSubCells
    (
        List<VtkSubCell>& cells,
        const label collapseDir
    ) const;
    void buildPrismCollapsedSubCells(List<VtkSubCell>& cells) const;

    static pointField orientPrismVertices(const pointField& rawVertices);

    static scalar prismOrientationMeasure
    (
        const UList<point>& localPoints,
        const FixedList<label, 8>& pointLabels
    );

    static void orientPrismSubCell
    (
        const UList<point>& localPoints,
        FixedList<label, 8>& pointLabels,
        const label nPoints
    );

public:

    dgEquiSpacedTessellation(const label n, const label pOrder);

    inline label n() const
    {
        return n_;
    }

    inline label pOrder() const
    {
        return pOrder_;
    }

    static dgCellType cellType(const label nPoints);

    CellInfo cellInfo
    (
        const cellShape& shape,
        const UList<point>& meshPoints
    ) const;

    const List<SamplePoint>& samplePoints(const dgCellType type) const;
    const List<SamplePoint>& collapsedSamplePoints
    (
        const dgCellType type,
        const label collapseDir
    ) const;

    void appendVtkCells
    (
        const dgCellType type,
        const UList<point>& localPoints,
        const label pointOffset,
        DynamicList<label>& connectivity,
        DynamicList<label>& offsets,
        DynamicList<unsigned char>& cellTypes
    ) const;

    void appendCollapsedVtkCells
    (
        const dgCellType type,
        const label collapseDir,
        const UList<point>& localPoints,
        const label pointOffset,
        DynamicList<label>& connectivity,
        DynamicList<label>& offsets,
        DynamicList<unsigned char>& cellTypes
    ) const;

    void appendVtkCells
    (
        const dgCellType type,
        const UList<point>& localPoints,
        const UList<label>& pointLabels,
        DynamicList<label>& connectivity,
        DynamicList<label>& offsets,
        DynamicList<unsigned char>& cellTypes
    ) const;

    void appendCollapsedVtkCells
    (
        const dgCellType type,
        const label collapseDir,
        const UList<point>& localPoints,
        const UList<label>& pointLabels,
        DynamicList<label>& connectivity,
        DynamicList<label>& offsets,
        DynamicList<unsigned char>& cellTypes
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
