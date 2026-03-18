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
    Mathematical utilities for DG Gauss fields.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGaussFieldMath_H
#define Foam_dgGaussFieldMath_H

#include "dgFaceGaussFieldLiteralScalarMathTmp.H"
#include "dgFaceGaussFieldScalarMathTmp.H"
#include "dgFaceGaussFieldVectorMathTmp.H"
#include "dgFaceGaussFieldTensorMathTmp.H"

#include "dgCellGaussFieldLiteralScalarMathTmp.H"
#include "dgCellGaussFieldScalarMathTmp.H"
#include "dgCellGaussFieldVectorMathTmp.H"
#include "dgCellGaussFieldTensorMathTmp.H"

#include "dgGaussFieldLiteralScalarMathTmp.H"
#include "dgGaussFieldScalarMathTmp.H"
#include "dgGaussFieldVectorMathTmp.H"
#include "dgGaussFieldTensorMathTmp.H"

#include "boundaryGaussFieldMath.H"

#endif

// ************************************************************************* //
