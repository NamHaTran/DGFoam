/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::boundaryGaussField

Description
    Lightweight boundary Gauss-point field container.

    Stores only a flat list of values on boundary-face Gauss points.
    The ordering is defined externally by faceGaussField::extractBCGhostState()
    and must be preserved when assigning values back with
    faceGaussField::assignBCGhostState().

SourceFiles
    boundaryGaussFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundaryGaussField_H
#define Foam_boundaryGaussField_H

#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class boundaryGaussField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::boundaryGaussField
 * \brief Flat container of Gauss-point values on boundary faces.
 *
 * \tparam Type Field value type stored at each boundary Gauss point.
 *
 * \details
 * boundaryGaussField stores a one-dimensional sequence of values defined on
 * boundary-face Gauss points. The mapping between list index and geometric
 * location is managed externally by the corresponding faceGaussField helpers.
 */
template<class Type>
class boundaryGaussField
{
private:

        /** \brief Stored values at boundary-face Gauss points. */
        List<Type> values_;


public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** \brief Construct an empty boundary Gauss field. */
        boundaryGaussField();

        /**
         * \brief Construct a boundary Gauss field with a given size.
         *
         * \param nSize Number of stored boundary Gauss-point values.
         */
        boundaryGaussField
        (
            const label nSize
        );

        /**
         * \brief Copy constructor.
         *
         * \param other Boundary Gauss field to copy from.
         */
        boundaryGaussField
        (
            const boundaryGaussField<Type>& other
        );


    //====================================================================//
    //                             Assignment
    //====================================================================//

        /**
         * \brief Copy assignment operator.
         *
         * \param other Boundary Gauss field to copy from.
         * \return Reference to this object.
         */
        boundaryGaussField<Type>& operator=
        (
            const boundaryGaussField<Type>& other
        );

        /**
         * \brief Assign a uniform value to all stored boundary Gauss points.
         *
         * \param value Uniform value assigned to every entry.
         * \return Reference to this object.
         */
        boundaryGaussField<Type>& operator=
        (
            const Type& value
        );


    //====================================================================//
    //                              Access
    //====================================================================//

        /** \brief Return the number of stored boundary Gauss-point values. */
        inline label size() const;


    //====================================================================//
    //                           Value access
    //====================================================================//

        /**
         * \brief Return writable access to one stored value.
         *
         * \param i Flat boundary Gauss-point index.
         * \return Reference to the requested value.
         */
        inline Type& operator[](const label i);

        /**
         * \brief Return read-only access to one stored value.
         *
         * \param i Flat boundary Gauss-point index.
         * \return Const reference to the requested value.
         */
        inline const Type& operator[](const label i) const;
};

} // End namespace Foam

#include "boundaryGaussFieldI.H"

#endif

// ************************************************************************* //
