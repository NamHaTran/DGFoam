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
    dgFaceType

Description
    Enumeration of supported DG face types.

    This enum is used in DG mesh handling to distinguish between
    triangular and quadrilateral faces in a topology-independent
    manner.

SourceFiles
    dgFaceType.H

\*---------------------------------------------------------------------------*/

#ifndef DGFACETYPE_H
#define DGFACETYPE_H

#include "IOstreams.H"

/**
 * \enum dgFaceType
 * \brief DG face topology identifier
 *
 * \details
 * Enumerates supported face types in the DG framework.
 * The value is typically associated with a DG face to
 * select appropriate basis functions and quadrature rules.
 */
enum class dgFaceType
{
    /** Triangular face */
    TRI,

    /** Quadrilateral face */
    QUAD,

    /** Invalid or unsupported face type */
    INVALID
};

/**
 * \brief Stream output operator for dgFaceType
 *
 * Writes a human-readable string representation of the face type.
 */
inline Foam::Ostream& operator<<(Foam::Ostream& os, const dgFaceType type)
{
    switch (type)
    {
        case dgFaceType::TRI:   return os << "TRI";
        case dgFaceType::QUAD:  return os << "QUAD";
        default:                return os << "UNKNOWN";
    }
}

#endif // DGFACETYPE_H
