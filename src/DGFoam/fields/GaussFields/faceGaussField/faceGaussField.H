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
// #include "faceGaussFieldVectorI.H"

// Forward declaration
class fieldsContext;

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
    const dgGeomMesh* mesh_;

    //- Cell ID
    label cellID_;

    //- Reference to the geometric cell
    const dgGeomCell* cell_;

    //- List of face IDs of the cell
    const labelList* facesID_;

    //- Number of faces
    label nFaces_;

    //- Number of Gauss points per face
    label nGaussPerFace_;

    //- Total number of Gauss points on all faces
    label nGauss_;

    //- List of pointers to dgGeomFace of this cell
    List<dgGeomFace*> faces_;

    //- List of DOF pointers for current cell and neighbor cells
    List<const cellDof<Type>*> cellsDof_;

    //- Per-face offset in Gauss values array
    List<label> gaussOffset_;

    //- Values on plus side (neighbor side) of faces
    List<Type> plusValues_;

    //- Values on minus side (owner side) of faces
    List<Type> minusValues_;

public:

    //- Default (empty) constructor
    faceGaussField();

    //- Construct from list of cell DOFs (P cell + neighbor cells) and mesh
    faceGaussField
    (
        List<const cellDof<Type>*>& cellsDof,
        const dgGeomMesh* mesh
    );

    //- Copy constructor
    faceGaussField(const faceGaussField<Type>& other);

    //- Construct from cellID, mesh, and initial values
    faceGaussField
    (
        const label cellID,
        const dgGeomMesh* mesh,
        const Type& initialValues
    );

    //- Construct from cellID and mesh (values will be filled later)
    faceGaussField
    (
        const label cellID,
        const dgGeomMesh* mesh
    );

    //- Set cell dof
    void setCellsDof(const List<const cellDof<Type>*>& cellsDof) {cellsDof_ = cellsDof;}

    //- Return number of faces
    inline label nFaces() const { return nFaces_; }

    //- Return number of Gauss points on all faces
    inline label nGauss() const { return nGauss_; }

    //- Interpolate values from DOF to Gauss points
    void interpolateFromDof();

    // Accessors

        //- Return cell ID
        inline label cellID() const { return cellID_; }

        //- Return pointer to geometric mesh
        inline const dgGeomMesh* dgMesh() const 
        {
            if (!mesh_)
            {
                FatalErrorInFunction
                    << "mesh_ pointer of faceGaussField is null" << abort(FatalError);
            }
            return mesh_;
        }

        //- Return value on plus side at global Gauss point index
        inline Type& plusValueAt(const label nG) { return plusValues_[nG]; }

        //- Return value on minus side at global Gauss point index
        inline Type& minusValueAt(const label nG) { return minusValues_[nG]; }

        //- Return value on plus side (const) at global Gauss point index
        inline const Type& plusValue(const label nG) const { return plusValues_[nG]; }

        //- Return value on minus side (const) at global Gauss point index
        inline const Type& minusValue(const label nG) const { return minusValues_[nG]; }

        //- Return value on plus side at local ID faceI and Gauss point index
        inline const Type& plusValue(const label faceI, const label gpI) const
        {
            return plusValues_[gaussOffset_[faceI] + gpI];
        }

        //- Return value on minus side at local ID faceI and Gauss point index
        inline const Type& minusValue(const label faceI, const label gpI) const
        {
            return minusValues_[gaussOffset_[faceI] + gpI];
        }

        //- Assign value on plus side at local ID faceI and Gauss point index
        inline Type& plusValueOnFace(const label faceI, const label gpI)
        {
            return plusValues_[gaussOffset_[faceI] + gpI];
        }

        //- Assign value on minus side at local ID faceI and Gauss point index
        inline Type& minusValueOnFace(const label faceI, const label gpI)
        {
            return minusValues_[gaussOffset_[faceI] + gpI];
        }

        // Operators

            //- Copy assignment operator
            faceGaussField<Type>& operator=(const faceGaussField<Type>& other);

            //- Assign a constant value to all Gauss points
            faceGaussField<Type>& operator=(const Type& value);

            //- Overload << operator for printing
            friend Ostream& operator<<(Ostream& os, const faceGaussField& fgf)
            {
                os << "Gauss field values for all faces of cell " << fgf.cellID() << nl;

                for (label i = 0; i < fgf.nGauss(); ++i)
                {
                    os << "  Gauss pt " << i << " on minus side: " << fgf.minusValue(i) << nl;
                    os << "  Gauss pt " << i << " on plus  side: " << fgf.plusValue(i) << nl;
                }

                return os;
            }

        // Specialisations for vector type
        faceGaussField<scalar> x() const;
        faceGaussField<scalar> y() const;
        faceGaussField<scalar> z() const;
        void setX(const faceGaussField<scalar>& fx);
        void setY(const faceGaussField<scalar>& fy);
        void setZ(const faceGaussField<scalar>& fz);
        void setXYZ
        (
            const faceGaussField<scalar>& fx,
            const faceGaussField<scalar>& fy,
            const faceGaussField<scalar>& fz
        );

        // Specialisations for scalar and vector type
        void limitMax(const scalar maxVal);
        void limitMin(const scalar minVal);
        void limitRange(const scalar minVal, const scalar maxVal);
};

} // End namespace Foam

// Specialisations for scalar type
#include "faceGaussFieldScalarI.H"

// Specialisations for vector type
#include "faceGaussFieldVectorI.H"

#endif

// ************************************************************************* //

