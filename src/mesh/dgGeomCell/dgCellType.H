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
    Foam::dgFacePosition

Description
    Enumeration for DG face positions on reference elements.

    The face positions are defined according to the vertex ordering
    convention used in:

        Sherwin, S.J., Karniadakis, G.E.
        "Spectral/hp Element Methods for Computational Fluid Dynamics",
        2nd Edition, Oxford University Press, 2005.

    The enumeration encodes the relative position of a face within a
    reference element using the same vertex labels (A–H) as in the
    literature, ensuring consistent evaluation of DG face integrals.

SourceFiles
    dgFacePosition.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFacePosition_H
#define Foam_dgFacePosition_H

#include "error.H"
#include "Ostream.H"

namespace Foam
{

/**
 * \enum dgFacePosition
 * \brief Face-position identifiers for a reference hexahedron
 *
 * \details
 * Identifies the position of a face on the reference hexahedral
 * element using counter-clockwise vertex ordering when viewed
 * from outside the element.
 */
enum class dgFacePosition
{
    /** Face ACEG */
    ACEG = 0,

    /** Face BDFH */
    BDFH,

    /** Face ABEF */
    ABEF,

    /** Face CDGH */
    CDGH,

    /** Face ABCD */
    ABCD,

    /** Face EFGH */
    EFGH,

    /** Number of valid face positions */
    NUM_POSITIONS,

    /** No face position (used for boundary faces) */
    NONE
};


/**
 * \enum dgFacePositionOnTet
 * \brief Face-position identifiers for a tetrahedral element
 *
 * \details
 * Face positions are defined by collapsing the reference hexahedron
 * onto a tetrahedron while preserving vertex ordering.
 */
enum class dgFacePositionOnTet
{
    /** Face BCE */
    BCE = 0,

    /** Face ACE */
    ACE,

    /** Face ABE */
    ABE,

    /** Face ABC */
    ABC,

    /** Number of faces */
    NUM_FACES,

    /** No face position (used for boundary faces) */
    NONE
};


/**
 * \enum dgFacePositionOnPyramid
 * \brief Face-position identifiers for a pyramidal element
 *
 * \details
 * Enumerates the base and triangular faces of the reference pyramid
 * element following DG vertex ordering conventions.
 */
enum class dgFacePositionOnPyramid
{
    /** Base face ABCD */
    ABCD = 0,

    /** Triangular face ACE */
    ACE,

    /** Triangular face CDE */
    CDE,

    /** Triangular face BDE */
    BDE,

    /** Triangular face ABE */
    ABE,

    /** Number of faces */
    NUM_FACES,

    /** No face position (used for boundary faces) */
    NONE
};


/**
 * \enum dgFacePositionOnPrism
 * \brief Face-position identifiers for a prismatic (wedge) element
 *
 * \details
 * Identifies triangular and quadrilateral faces of the reference
 * prism element using DG-consistent vertex ordering.
 */
enum class dgFacePositionOnPrism
{
    /** Triangular face ABE */
    ABE = 0,

    /** Triangular face CDF */
    CDF,

    /** Quadrilateral face BCEF */
    BCEF,

    /** Quadrilateral face ABCD */
    ABCD,

    /** Quadrilateral face ADEF */
    ADEF,

    /** Number of faces */
    NUM_FACES,

    /** No face position (used for boundary faces) */
    NONE
};


/**
 * \brief Map tetrahedron face position to reference hexahedron face position
 */
inline dgFacePosition mapFacePositionFromTet(const dgFacePositionOnTet f)
{
    switch (f)
    {
        case dgFacePositionOnTet::ACE:  return dgFacePosition::ACEG;
        case dgFacePositionOnTet::BCE:  return dgFacePosition::BDFH;
        case dgFacePositionOnTet::ABE:  return dgFacePosition::ABEF;
        case dgFacePositionOnTet::ABC:  return dgFacePosition::ABCD;
        case dgFacePositionOnTet::NONE: return dgFacePosition::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePositionOnTet enum"
                                 << abort(FatalError);
    }

    return dgFacePosition::ABCD;
}


/**
 * \brief Map reference hexahedron face position to tetrahedron face position
 */
inline dgFacePositionOnTet mapFacePositionToTet(const dgFacePosition f)
{
    switch (f)
    {
        case dgFacePosition::ACEG: return dgFacePositionOnTet::ACE;
        case dgFacePosition::BDFH: return dgFacePositionOnTet::BCE;
        case dgFacePosition::ABEF: return dgFacePositionOnTet::ABE;
        case dgFacePosition::ABCD: return dgFacePositionOnTet::ABC;
        case dgFacePosition::NONE: return dgFacePositionOnTet::NONE;
        default:
            FatalErrorInFunction << "Invalid dgFacePosition enum"
                                 << abort(FatalError);
    }

    return dgFacePositionOnTet::ABC;
}


/**
 * \brief Stream output operator for dgFacePosition
 */
inline Ostream& operator<<(Ostream& os, const dgFacePosition type)
{
    switch (type)
    {
        case dgFacePosition::ACEG: return os << "ACEG";
        case dgFacePosition::BDFH: return os << "BDFH";
        case dgFacePosition::ABEF: return os << "ABEF";
        case dgFacePosition::CDGH: return os << "CDGH";
        case dgFacePosition::ABCD: return os << "ABCD";
        case dgFacePosition::EFGH: return os << "EFGH";
        case dgFacePosition::NONE: return os << "NONE";
        default:                   return os << "UNKNOWN";
    }
}

} // End namespace Foam

#endif // Foam_dgFacePosition_H

// ************************************************************************* //
