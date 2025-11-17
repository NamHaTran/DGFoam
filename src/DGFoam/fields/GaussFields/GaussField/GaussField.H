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

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class GaussField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GaussField
{
private:

    //- Pointer to DoF field (read-only, non-owning)
    const dofField<Type>* dofField_;

    //- Cell ID
    label cellID_;

    //- Pointer to DG mesh (non-owning)
    const dgGeomMesh* mesh_;

    //- Cell Gauss values
    cellGaussField<Type> cellField_;

    //- Face Gauss values
    faceGaussField<Type> faceField_;

public:

    // Constructors ----------------------------------------------------------

    //- Default constructor
    GaussField();

    //- Construct from dofField, cellID and mesh
    GaussField
    (
        const dofField<Type>* dofField,
        label cellID,
        const dgGeomMesh* mesh
    );

    //- Construct container only
    GaussField(label cellID, const dgGeomMesh* mesh);

    //- Construct container with initial value
    GaussField(label cellID, const dgGeomMesh* mesh, const Type& initVal);

    //- Copy constructor
    GaussField(const GaussField<Type>& other);


    // Assignment ------------------------------------------------------------

    GaussField<Type>& operator=(const GaussField<Type>& other);
    GaussField<Type>& operator=(const Type& value);


    // Access ---------------------------------------------------------------

    inline const dgGeomMesh* dgMesh() const;

    inline label cellID() const { return cellID_; }

    inline const dofField<Type>* dofFieldPtr() const { return dofField_; }

    inline cellGaussField<Type>& cellField()       { return cellField_; }
    inline const cellGaussField<Type>& cellField() const { return cellField_; }

    inline faceGaussField<Type>& faceField()       { return faceField_; }
    inline const faceGaussField<Type>& faceField() const { return faceField_; }


    // High-level operations ------------------------------------------------

    //- Interpolate values from DOF into Gauss points
    void interpolateFromDof();


    // Direct Gauss access --------------------------------------------------

    inline Type& cellValueAt(label gpI);
    inline const Type& cellValueAt(label gpI) const;

    inline Type& facePlusValueAt(label faceI, label gpI);
    inline Type& faceMinusValueAt(label faceI, label gpI);

    inline const Type& facePlusValueAt(label faceI, label gpI) const;
    inline const Type& faceMinusValueAt(label faceI, label gpI) const;


    // Pretty print
    friend Ostream& operator<<(Ostream& os, const GaussField& gf)
    {
        os << "GaussField for cell " << gf.cellID() << ":\n";
        os << "Cell Gauss Field:\n" << gf.cellField() << "\n";
        os << "Face Gauss Field:\n" << gf.faceField() << "\n";
        return os;
    }


    // Vector specialisations
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

    // Limiters
    void limitMax(scalar maxVal);
    void limitMin(scalar minVal);
    void limitRange(scalar minVal, scalar maxVal);
};

} // End namespace Foam


// * * * * * * * * Template implementation * * * * * * * *
#include "GaussFieldI.H"
#include "GaussFieldVectorI.H"
#include "GaussFieldScalarI.H"

#endif


