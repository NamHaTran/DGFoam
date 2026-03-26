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

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundaryGaussField_H
#define Foam_boundaryGaussField_H

#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class boundaryGaussField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class boundaryGaussField
{
private:

        //- Stored Gauss-point values on boundary faces
        List<Type> values_;


public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        //- Default constructor
        boundaryGaussField();

        //- Construct boundary field with a given size
        boundaryGaussField
        (
            const label nSize
        );

        //- Copy constructor
        boundaryGaussField
        (
            const boundaryGaussField<Type>& other
        );


    //====================================================================//
    //                             Assignment
    //====================================================================//

        //- Copy assignment
        boundaryGaussField<Type>& operator=
        (
            const boundaryGaussField<Type>& other
        );

        //- Assign uniform value
        boundaryGaussField<Type>& operator=
        (
            const Type& value
        );


    //====================================================================//
    //                              Access
    //====================================================================//

        //- Return number of stored values
        inline label size() const;


    //====================================================================//
    //                           Value access
    //====================================================================//

        //- Mutable element access
        inline Type& operator[](const label i);

        //- Const element access
        inline const Type& operator[](const label i) const;
};

} // End namespace Foam

#include "boundaryGaussFieldI.H"

#endif

// ************************************************************************* //
