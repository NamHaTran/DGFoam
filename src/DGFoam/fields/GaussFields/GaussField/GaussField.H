/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024-2025 Ha Nam Tran
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
    Foam::faceGaussField

Description
    Wrapper class to store and manage field values at Gauss integration points
    within a single cell and its faces.

    This class provides access to both internal cell Gauss values and the Gauss
    values located on each face of the cell. It is designed to be constructed
    either from a pointer to a DoF field (e.g. dofField<Type>) or simply as a
    container to hold computed values.

    - The cell interior values are stored in a cellGaussField<Type> instance.
    - The face values (both plus and minus sides) are stored in a faceGaussField<Type>.
    - Supports various initialization modes: from DoFs, from constant value, or empty.

SourceFiles
    faceGaussField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GaussField_H
#define Foam_GaussField_H

#include "cellGaussField.H"
#include "faceGaussField.H"
#include "dofField.H"
#include "dgGeomMesh.H"
// #include "GaussFieldVectorI.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class GaussField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GaussField
{
private:

    //- Pointer to the DoF field (read-only)
    const dofField<Type>* dofField_;

    //- Cell ID associated with this field
    label cellID_;

    //- Reference to the DG mesh
    const dgGeomMesh* mesh_;

    //- Gauss field values inside the cell
    cellGaussField<Type> cellField_;

    //- Gauss field values on faces of the cell
    faceGaussField<Type> faceField_;

public:

    // Constructors

        //- Default (empty) constructor
        GaussField();

        //- Construct from dofField, cellID and mesh
        GaussField(const dofField<Type>* dofField, label cellID, const dgGeomMesh* mesh);

        //- Construct as container only (empty data)
        GaussField(label cellID, const dgGeomMesh* mesh);

        //- Construct as container with initial value
        GaussField(label cellID, const dgGeomMesh* mesh, const Type& initVal);

        //- Copy constructor
        GaussField(const GaussField<Type>& other);

    //- Assignment operator
    GaussField<Type>& operator=(const GaussField<Type>& other);

    //- Assign a constant value to all Gauss points
    GaussField<Type>& operator=(const Type& value);

    // Member access

        //- Return const pointer to geometric mesh
        inline const dgGeomMesh* dgMesh() const 
        {
            if (!mesh_)
            {
                FatalErrorInFunction
                    << "mesh_ pointer of GaussField is null" << abort(FatalError);
            }
            return mesh_;
        }

        //- Return reference to cellGaussField
        inline cellGaussField<Type>& cellField() { return cellField_; }

        //- Return const reference to cellGaussField
        inline const cellGaussField<Type>& cellField() const { return cellField_; }

        //- Return reference to faceGaussField
        inline faceGaussField<Type>& faceField() { return faceField_; }

        //- Return const reference to faceGaussField
        inline const faceGaussField<Type>& faceField() const { return faceField_; }

        //- Return cell ID
        inline label cellID() const { return cellID_; }

        //- Return const pointer to associated dofField (may be nullptr)
        inline const dofField<Type>* dofFieldPtr() const { return dofField_; }

    // Member Functions

        //- Interpolate values at Gauss points from DOF using basis functions
        void interpolateFromDof();

        //- Overload << operator for printing
        friend Ostream& operator<<(Ostream& os, const GaussField& gf)
        {
            os << "GaussField for cell " << gf.cellID() << ":\n";
            os << "Cell Gauss Field:\n" << gf.cellField() << "\n";
            os << "Face Gauss Field:\n" << gf.faceField() << "\n";
            return os;
        }

        //- Assign value to volume Gauss point
        inline Type& cellValueAt(const label nG)
        {
            return cellField_[nG];
        } 

        //- Access value at volume Gauss point (const)
        inline const Type& cellValueAt(const label nG) const
        {
            return cellField_[nG];
        } 

        //- Assign value to face Gauss point (plus side), faceI is local face ID on list faces of the cell
        inline Type& facePlusValueAt(const label faceI, const label gpI)
        {
            return faceField_.plusValueOnFace(faceI, gpI);
        }

        //- Access value at face Gauss point (plus side, const), faceI is local face ID on list faces of the cell
        inline const Type& facePlusValueAt(const label faceI, const label gpI) const
        {
            return faceField_.plusValue(faceI, gpI);
        }

        //- Assign value to face Gauss point (minus side), faceI is local face ID on list faces of the cell
        inline Type& faceMinusValueAt(const label faceI, const label gpI)
        {
            return faceField_.minusValueOnFace(faceI, gpI);
        }

        //- Access value at face Gauss point (minus side, const), faceI is local face ID on list faces of the cell
        inline const Type& faceMinusValueAt(const label faceI, const label gpI) const
        {
            return faceField_.minusValue(faceI, gpI);
        }

        //- Return cellGaussField
        inline cellGaussField<Type>& getCellGaussField() { return cellField_; }

        //- Return faceGaussField
        inline faceGaussField<Type>& getFaceGaussField() { return faceField_; }

    // Specialisations for vector type
    GaussField<scalar> x() const;
    GaussField<scalar> y() const;
    GaussField<scalar> z() const;
    void setX(const GaussField<scalar>& fx);
    void setY(const GaussField<scalar>& fy);
    void setZ(const GaussField<scalar>& fz);
    void setXYZ
    (
        const GaussField<scalar>& fx,
        const GaussField<scalar>& fy,
        const GaussField<scalar>& fz
    );

    // Specialisations for scalar and vector type
    void limitMax(const scalar maxVal);
    void limitMin(const scalar minVal);
    void limitRange(const scalar minVal, const scalar maxVal);
};

} // End namespace Foam

// Specialisations for vector type
#include "GaussFieldVectorI.H"

#endif

// ************************************************************************* //

