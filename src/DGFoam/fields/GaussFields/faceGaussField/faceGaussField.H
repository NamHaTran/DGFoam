/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR by YOUR NAME
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
    Stores Gauss point values on all faces of a cell for DG method.
    Supports scalar, vector, and tensor values.
    
    - Stores both minus and plus side values of all Gauss points.
    - Designed for cell-centric DG formulation.
    - Uses DOF from the current cell and all neighboring cells.

SourceFiles
    faceGaussField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceGaussField_H
#define Foam_faceGaussField_H

#include "dgGeomMesh.H"
#include "dgGeomCell.H"
#include "dgGeomFace.H"
#include "cellDof.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceGaussField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faceGaussField
{
private:

    //- Reference to geometric mesh
    const dgGeomMesh& mesh_;

    //- Cell ID
    label cellID_;

    //- Pointer to the geometric cell
    dgGeomCell* cell_;

    //- List of face IDs of the cell
    const labelList& facesID_;

    //- Number of faces
    label nFaces_;

    //- Number of Gauss points per face
    label nGaussPerFace_;

    //- Total number of Gauss points on all faces
    label nGauss_;

    //- List of pointers to dgGeomFace of this cell
    List<dgGeomFace*> faces_;

    //- List of DOF pointers for current cell and neighbor cells
    List<cellDof<Type>*> cellsDof_;

    //- Per-face offset in Gauss values array
    List<label> gaussOffset_;

    //- Values on plus side (neighbor side) of faces
    List<Type> plusValues_;

    //- Values on minus side (owner side) of faces
    List<Type> minusValues_;

public:

    //- Construct from list of cell DOFs (P cell + neighbor cells) and mesh
    faceGaussField
    (
        const List<cellDof<Type>*>& cellsDof,
        const dgGeomMesh& mesh
    );

    //- Construct from list of cell DOFs, mesh, and OpenFOAM volField
    faceGaussField
    (
        const List<cellDof<Type>*>& cellsDof,
        const GeometricField<Type, fvPatchField, volMesh>& foamField,
        const dgGeomMesh& mesh
    );

    //- Construct from list of cell DOFs, mesh, and initial values
    faceGaussField
    (
        const List<cellDof<Type>*>& cellsDof,
        const dgGeomMesh& mesh,
        const Type& initialValues
    );
    
    //- Copy constructor
    faceGaussField(const faceGaussField<Type>& other);

    //- Construct from cellID, mesh, and initial values
    faceGaussField
    (
        const label cellID,
        const dgGeomMesh& mesh,
        const Type& initialValues
    );

    //- Construct from cell pointer, mesh
    // This is empty and will be filled later
    faceGaussField
    (
        const label cellID,
        const dgGeomMesh& mesh
    );

    //- Return number of faces
    inline label nFaces() const { return nFaces_; }

    //- Return number of Gauss points on all faces
    inline label nGauss() const { return nGauss_; }

    //- Interpolate values from DOF to Gauss points
    void interpolate();

    //- Accessors:

        //- Return cell ID
        inline label cellID() const {return cellID_;};

        //- Return const reference to geometric mesh (for const object)
        inline const dgGeomMesh& dgMesh() const {return mesh_;};

        //- Return value on plus side at global Gauss point index, allow modifying value
        inline Type& plusValueAt(const label nG) {return plusValues_[nG];};

        //- Return value on minus side at global Gauss point index, allow modifying value
        inline Type& minusValueAt(const label nG) {return minusValues_[nG];};

        //- Return value on plus side by global Gauss point ID
        inline const Type& plusValue(const label nG) const {return plusValues_[nG];};

        //- Return value on minus side by global Gauss point ID
        inline const Type& minusValue(const label nG) const {return minusValues_[nG];};

        //- Return value on plus side at faceI and Gauss point index
        inline const Type& plusValue(const label faceI, const label gpI) const {return plusValues_[gaussOffset_[faceI] + gpI];};

        //- Return value on minus side at faceI and Gauss point index
        inline const Type& minusValue(const label faceI, const label gpI) const {return minusValues_[gaussOffset_[faceI] + gpI];};

};

} // End namespace Foam

#ifdef NoRepository
    #include "faceGaussField.C"
#endif

#endif

// ************************************************************************* //
