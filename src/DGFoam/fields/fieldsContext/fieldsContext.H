/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C)
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
    Foam::fieldsContext

Description
    Context class to register and lookup cellGaussField and faceGaussField
    using string names. Stores internal pointers so fields can be shared and
    accessed during runtime with updated values.

    Supports multiple data types (scalar, vector, tensor, etc.) via template
    mechanism.

    Example usage:

        fieldsContext ctx;

        cellGaussField<vector> U(...);
        faceGaussField<scalar> TWall(...);

        ctx.regis(U, "U");
        ctx.regis(TWall, "TWall");

        const cellGaussField<vector>& Uref = ctx.lookupCellField<vector>("U");
        const faceGaussField<scalar>& TwRef = ctx.lookupFaceField<scalar>("TWall");

SourceFiles
    fieldsContext.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fieldsContext_H
#define Foam_fieldsContext_H

#include "cellGaussField.H"
#include "faceGaussField.H"
#include "HashTable.H"
#include "word.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fieldsContext Declaration
\*---------------------------------------------------------------------------*/

class fieldsContext
{
public:

    //- Field location: cell or face
    enum fieldType
    {
        CELL_FIELD,
        FACE_FIELD
    };

private:

    //- Template alias for cell field map
    template<class Type>
    using CellFieldMap = HashTable<const cellGaussField<Type>*>;

    //- Template alias for face field map
    template<class Type>
    using FaceFieldMap = HashTable<const faceGaussField<Type>*>;

    //- Generic maps storing pointer-to-maps using type string as key
    HashTable<void*> cellFieldMaps_;
    HashTable<void*> faceFieldMaps_;

    //- Helper: generate key string from Type
    template<class Type>
    static word typeKey();

public:

    //- Default constructor
    fieldsContext() = default;

    //- Register a cellGaussField
    template<class Type>
    void regis(const cellGaussField<Type>& field, const word& name);

    //- Register a faceGaussField
    template<class Type>
    void regis(const faceGaussField<Type>& field, const word& name);

    //- Lookup cell field
    template<class Type>
    const cellGaussField<Type>& lookupCellField(const word& name) const;

    //- Lookup face field
    template<class Type>
    const faceGaussField<Type>& lookupFaceField(const word& name) const;
};

} // End namespace Foam

#ifdef NoRepository
    #include "fieldsContext.C"
#endif

#endif

// ************************************************************************* //
