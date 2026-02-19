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

Class

Description
    Just a cell type enumeration for use in DG mesh handling.

See also

SourceFiles
    dgCellType.H

\*---------------------------------------------------------------------------*/

#ifndef DGCELLTYPE_H
#define DGCELLTYPE_H

#include "IOstreams.H"

enum class dgCellType
{
    TET,
    HEX,
    PRISM,
    PYRAMID,
    NONE, // No cell type (used for boundary faces)
    INVALID
};

//- Overwrite << operator
inline Foam::Ostream& operator<<(Foam::Ostream& os, const dgCellType type)
{
    switch (type)
    {
        case dgCellType::TET:   return os << "TET";
        case dgCellType::HEX:   return os << "HEX";
        case dgCellType::PRISM: return os << "PRISM";
        case dgCellType::PYRAMID: return os << "PYRAMID";
        case dgCellType::NONE:  return os << "NONE";
        default:                return os << "UNKNOWN";
    }
}
#endif // DGCELLTYPE_H
