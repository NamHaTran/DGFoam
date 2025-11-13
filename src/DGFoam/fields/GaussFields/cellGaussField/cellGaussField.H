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
    const dgGeomMesh* dgMesh_ = nullptr;

    //- Reference to the geometric cell
    const dgGeomCell* cell_ = nullptr;

    //- Number of Gauss points in this cell
    label nGauss_;

    //- Number of DOFs in this cell
    label nDof_;

    //- DOF values in this cell
    const cellDof<Type>* dof_;

    //- Field values at Gauss points
    List<Type> values_;

public:

    //- Default (empty) constructor
    cellGaussField();

    //- Constructor from DOF object
    cellGaussField
    (
        const dgGeomMesh* dgMesh,
        const cellDof<Type>* dof
    );

    //- Copy constructor
    cellGaussField(const cellGaussField<Type>& other);

    //- Initial value constructor
    cellGaussField
    (
        const label cellID,
        const dgGeomMesh* dgMesh,
        const Type& initVal
    );

    //- Constructor without initial values
    cellGaussField
    (
        const label cellID,
        const dgGeomMesh* dgMesh
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
    inline const dgGeomMesh* dgMesh() const 
    {
        if (!dgMesh_)
        {
            FatalErrorInFunction
                << "dgMesh_ pointer of cellGaussField is null" << abort(FatalError);
        }
        return dgMesh_;
    }

    //- Interpolate values at Gauss points from DOF using basis functions
    void interpolateFromDof();

    //- Set the geometric mesh pointer
    inline void setMesh(const dgGeomMesh& mesh)
    {
        if (dgMesh_)
        {
            WarningInFunction
                << "dgMesh_ pointer is already set. Overwriting existing pointer." 
                << nl << endl;
        }

        dgMesh_ = &mesh;
    }

    // Operators

        //- Copy assignment operator
        cellGaussField<Type>& operator=(const cellGaussField<Type>& other);

        //- Assign a constant value to all Gauss points
        cellGaussField<Type>& operator=(const Type& value);

        //- Overload << operator for printing
        friend Ostream& operator<<(Ostream& os, const cellGaussField& cgf)
        {
            os << "cellGaussField for cell " << cgf.cellID() << ":\n";
            for (label i = 0; i < cgf.size(); ++i)
            {
                os << "  Gauss point " << i << ": " << cgf[i] << "\n";
            }
            return os;
        }

    // Specialisations for vector type
    cellGaussField<scalar> x() const;
    cellGaussField<scalar> y() const;
    cellGaussField<scalar> z() const;
    void setX(const cellGaussField<scalar>& fx);
    void setY(const cellGaussField<scalar>& fy);
    void setZ(const cellGaussField<scalar>& fz);
    void setXYZ
    (
        const cellGaussField<scalar>& fx,
        const cellGaussField<scalar>& fy,
        const cellGaussField<scalar>& fz
    );

    // Specialisations for scalar and vector type
    void limitMax(const scalar maxVal);
    void limitMin(const scalar minVal);
    void limitRange(const scalar minVal, const scalar maxVal);
};

} // End namespace Foam

// Specialisations for scalar type
#include "cellGaussFieldScalarI.H"

// Specialisations for vector type
#include "cellGaussFieldVectorI.H"

#endif

// ************************************************************************* //

