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

\*---------------------------------------------------------------------------*/

#include "fieldsContext.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * * Type key generation * * * * * * * * * * * * * //
// Generate a unique key for each type used in fieldsContext
template<> word fieldsContext::typeKey<scalar>() { return "scalar"; }
template<> word fieldsContext::typeKey<vector>() { return "vector"; }
template<> word fieldsContext::typeKey<tensor>() { return "tensor"; }
template<> word fieldsContext::typeKey<symmTensor>() { return "symmTensor"; }
template<> word fieldsContext::typeKey<sphericalTensor>() { return "sphericalTensor"; }


// * * * * * * * * * * * * * * Register cell field * * * * * * * * * * * * * //

template<class Type>
void fieldsContext::regis(const cellGaussField<Type>& field, const word& name)
{
    word key = typeKey<Type>();

    CellFieldMap<Type>* mapPtr = nullptr;

    if (!cellFieldMaps_.found(key))
    {
        mapPtr = new CellFieldMap<Type>();
        cellFieldMaps_.insert(key, mapPtr);
    }
    else
    {
        mapPtr = static_cast<CellFieldMap<Type>*>(cellFieldMaps_[key]);
    }

    mapPtr->insert(name, &field);
}


// * * * * * * * * * * * * * * Register face field * * * * * * * * * * * * * //

template<class Type>
void fieldsContext::regis(const faceGaussField<Type>& field, const word& name)
{
    word key = typeKey<Type>();

    FaceFieldMap<Type>* mapPtr = nullptr;

    if (!faceFieldMaps_.found(key))
    {
        mapPtr = new FaceFieldMap<Type>();
        faceFieldMaps_.insert(key, mapPtr);
    }
    else
    {
        mapPtr = static_cast<FaceFieldMap<Type>*>(faceFieldMaps_[key]);
    }

    mapPtr->insert(name, &field);
}


// * * * * * * * * * * * * * * Lookup cell field * * * * * * * * * * * * * * //

template<class Type>
const cellGaussField<Type>& fieldsContext::lookupCellField(const word& name) const
{
    word key = typeKey<Type>();

    auto iter = cellFieldMaps_.find(key);
    if (iter == cellFieldMaps_.end())
    {
        FatalErrorInFunction << "No cell field map for type: " << key << exit(FatalError);
    }

    const CellFieldMap<Type>* mapPtr = static_cast<const CellFieldMap<Type>*>(iter());

    const auto fIter = mapPtr->find(name);
    if (fIter == mapPtr->end())
    {
        FatalErrorInFunction << "No cellGaussField<" << key << "> named '" << name << "' found" << exit(FatalError);
    }

    return *fIter();
}


// * * * * * * * * * * * * * * Lookup face field * * * * * * * * * * * * * * //

template<class Type>
const faceGaussField<Type>& fieldsContext::lookupFaceField(const word& name) const
{
    word key = typeKey<Type>();

    auto iter = faceFieldMaps_.find(key);
    if (iter == faceFieldMaps_.end())
    {
        FatalErrorInFunction << "No face field map for type: " << key << exit(FatalError);
    }

    const FaceFieldMap<Type>* mapPtr = static_cast<const FaceFieldMap<Type>*>(iter());

    const auto fIter = mapPtr->find(name);
    if (fIter == mapPtr->end())
    {
        FatalErrorInFunction << "No faceGaussField<" << key << "> named '" << name << "' found" << exit(FatalError);
    }

    return *fIter();
}


// * * * * * * * * * * * * * * Template instantiation * * * * * * * * * * * * //

// Instantiate the template functions for all types used in fieldsContext
// This is necessary to ensure that the template functions are compiled and linked correctly.
// Note: This is a common practice in C++ to avoid linker errors for template functions.
// If you add new types, you should also add them here.
// If you remove a type, you should also remove it from here.
#define makeFieldsContextFuncs(Type) \
    template void fieldsContext::regis<Type>(const cellGaussField<Type>&, const word&); \
    template void fieldsContext::regis<Type>(const faceGaussField<Type>&, const word&); \
    template const cellGaussField<Type>& fieldsContext::lookupCellField<Type>(const word&) const; \
    template const faceGaussField<Type>& fieldsContext::lookupFaceField<Type>(const word&) const;

// Register all types used in fieldsContext using the macro
makeFieldsContextFuncs(scalar)
makeFieldsContextFuncs(vector)
makeFieldsContextFuncs(tensor)
makeFieldsContextFuncs(symmTensor)
makeFieldsContextFuncs(sphericalTensor)

// Delete the macro to avoid redefinition errors
#undef makeFieldsContextFuncs

} // End namespace Foam

// ************************************************************************* //
