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

    //- Pointer to geometric mesh (non-owning)
    const dgGeomMesh* mesh_;

    //- Cell ID associated with this face field
    label cellID_;

    //- Pointer to geometric cell (non-owning)
    const dgGeomCell* cell_;

    //- Pointer to face ID list (non-owning)
    const labelList* facesID_;

    //- Number of faces in the cell
    label nFaces_;

    //- Number of Gauss points per face
    label nGaussPerFace_;

    //- Total number of Gauss points on all faces
    label nGauss_;

    //- List of face pointers (non-owning)
    List<dgGeomFace*> faces_;

    //- Pointer list of DOFs for owner and neighbour cells
    List<const cellDof<Type>*> cellsDof_;

    //- Per-face Gauss offset index
    List<label> gaussOffset_;

    //- Per-face patchID (if any)
    List<label> ownerPatchIDs_;

    //- Values on plus side (neighbour)
    List<Type> plusValues_;

    //- Values on minus side (owner)
    List<Type> minusValues_;

public:

    // Constructors ------------------------------------------------------------

    //- Default constructor
    faceGaussField();

    //- Construct from P-cell and neighbour DOF pointers
    faceGaussField
    (
        List<const cellDof<Type>*>& cellsDof,
        const dgGeomMesh* mesh
    );

    //- Copy constructor
    faceGaussField(const faceGaussField<Type>& other);

    //- Construct and set uniform initial values
    faceGaussField
    (
        const label cellID,
        const dgGeomMesh* mesh,
        const Type& initVal
    );

    //- Construct empty, fill values later
    faceGaussField
    (
        const label cellID,
        const dgGeomMesh* mesh
    );

    // Access -----------------------------------------------------------------

    inline label cellID() const { return cellID_; }

    inline const dgGeomMesh* dgMesh() const
    {
        if (!mesh_)
        {
            FatalErrorInFunction << "mesh_ pointer is null" << abort(FatalError);
        }
        return mesh_;
    }

    inline label nFaces() const { return nFaces_; }
    inline label nGauss() const { return nGauss_; }
    inline label nGaussPerFace() const { return nGaussPerFace_; }

    // Access to list of dgGeomFace
    inline const List<dgGeomFace*>& faces() const
    {
        return faces_;
    }

    // Access to value on plus side using its global Gauss id
    inline Type& plusValueAt(const label gi)       { return plusValues_[gi]; }

    // Access to value on minus side using its global Gauss id
    inline Type& minusValueAt(const label gi)      { return minusValues_[gi]; }

    // Access to value on plus side using its global Gauss id (const)
    inline const Type& plusValue(const label gi) const  { return plusValues_[gi]; }

    // Access to value on minus side using its global Gauss id (const)
    inline const Type& minusValue(const label gi) const { return minusValues_[gi]; }

    // Access to value of plus side base on local face id and local Gauss point ID on that face
    inline Type& plusValueOnFace(const label fI, const label gI)
    {
        return plusValues_[gaussOffset_[fI] + gI];
    }

    inline const Type& plusValueOnFace(const label fI, const label gI) const
    {
        return plusValues_[gaussOffset_[fI] + gI];
    }

    // Access to value of minus side base on local face id and local Gauss point ID on that face
    inline Type& minusValueOnFace(const label fI, const label gI)
    {
        return minusValues_[gaussOffset_[fI] + gI];
    }

    inline const Type& minusValueOnFace(const label fI, const label gI) const
    {
        return minusValues_[gaussOffset_[fI] + gI];
    }

    // Check whether the face is boundary
    inline bool isBoundary(const label fI) const
    {
        if (ownerPatchIDs_[fI] >= 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    // Get owner patch's id of face
    inline label getOwnerPatchID(const label fI) const
    {
        return ownerPatchIDs_[fI];
    }

    // Set DOF pointers
    inline void setCellsDof(const List<const cellDof<Type>*>& lst)
    {
        cellsDof_ = lst;
    }

    // Member Functions --------------------------------------------------------

    //- Interpolate DOF → Gauss values
    void interpolateFromDof();

    // Operators ---------------------------------------------------------------

    //- Assignment operator
    faceGaussField<Type>& operator=(const faceGaussField<Type>& other);

    //- Assign uniform value
    faceGaussField<Type>& operator=(const Type& v);

    // Printing ---------------------------------------------------------------

    friend Ostream& operator<<
    (
        Ostream& os,
        const faceGaussField<Type>& fgf
    )
    {
        os << "faceGaussField for cell " << fgf.cellID() << nl;

        for (label i = 0; i < fgf.nGauss(); ++i)
        {
            os << "  Gauss pt " << i << " minus: " << fgf.minusValue(i) << nl;
            os << "  Gauss pt " << i << " plus : " << fgf.plusValue(i)  << nl;
        }
        return os;
    }

    // Specialisation of vector/scalar components
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

    // Limiting
    void limitMax(const scalar maxVal);
    void limitMin(const scalar minVal);
    void limitRange(const scalar minVal, const scalar maxVal);
};

} // End namespace Foam

#include "faceGaussFieldI.H"
#include "faceGaussFieldScalarI.H"
#include "faceGaussFieldVectorI.H"

#endif

// ************************************************************************* //


