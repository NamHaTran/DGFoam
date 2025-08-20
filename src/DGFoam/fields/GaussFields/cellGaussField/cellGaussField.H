/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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
    Foam::cellGaussField

Description
    Field values at Gauss points in a DG cell.

SourceFiles
    cellGaussField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellGaussField_H
#define Foam_cellGaussField_H

#include "label.H"
#include "List.H"
#include "dgGeomMesh.H"
#include "dgGeomCell.H"
#include "cellDof.H"
#include "fieldsContext.H"

// Forward declaration
class fieldsContext;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellGaussField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellGaussField
{
private:

    //- Cell ID in the mesh
    label cellID_;

    //- Reference to geometric mesh
    const dgGeomMesh& dgMesh_;

    //- Reference to the geometric cell
    const dgGeomCell& cell_;

    //- Number of Gauss points in this cell
    label nGauss_;

    //- Number of DOFs in this cell
    label nDof_;

    //- DOF values in this cell
    List<Type> dof_;

    //- Field values at Gauss points
    List<Type> values_;

    //- Pointer to the fields context for registration
    const fieldsContext* ctxPtr_ = nullptr;

public:

    //- Constructor from cell ID and mesh
    cellGaussField
    (
        const label cellID, 
        const dgGeomMesh& dgMesh,
        const List<Type>& values
    );

    //- Constructor from DOF object
    cellGaussField(const cellDof<Type>& dof);

    //- Copy constructor
    cellGaussField(const cellGaussField<Type>& other);

    //- Initial value constructor
    cellGaussField
    (
        const label cellID,
        const dgGeomMesh& dgMesh,
        const Type& initVal
    );

    //- Constructor without initial values
    cellGaussField
    (
        const label cellID,
        const dgGeomMesh& dgMesh
    );

    //- Return reference to value at Gauss point
    inline Type& operator[](const label gpI) { return values_[gpI]; }

    //- Const access to value at Gauss point
    inline const Type& operator[](const label gpI) const { return values_[gpI]; }

    //- Return number of Gauss points
    inline label size() const { return nGauss_; }

    //- Return cell ID
    inline label cellID() const { return cellID_; }

    //- Return const reference to geometric mesh
    inline const dgGeomMesh& dgMesh() const { return dgMesh_; }

    //- Interpolate values at Gauss points from DOF using basis functions
    void interpolateFromDof();

    //- Set the fields context
    inline void setContext(const fieldsContext& ctx) { ctxPtr_ = &ctx;}

    //- Debug print
    void print() const;

    // Operators

        //- Copy assignment operator
        cellGaussField<Type>& operator=(const cellGaussField<Type>& other);

        //- Assign a constant value to all Gauss points
        cellGaussField<Type>& operator=(const Type& value);
};

} // End namespace Foam

#endif

// ************************************************************************* //

