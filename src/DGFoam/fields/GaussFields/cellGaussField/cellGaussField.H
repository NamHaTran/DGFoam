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

    //- Pointer to geometric mesh (non-owning)
    const dgGeomMesh* dgMesh_;

    //- Pointer to geometric cell (non-owning)
    const dgGeomCell* cell_;

    //- Number of Gauss points in this cell
    label nGauss_;

    //- Number of DOFs in this cell
    label nDof_;

    //- Pointer to DOF container (non-owning)
    const cellDof<Type>* dof_;

    //- Values at Gauss points
    List<Type> values_;


public:

    // Constructors ------------------------------------------------------------

    //- Default constructor
    cellGaussField();

    //- Construct from DG mesh + DOF pointer
    cellGaussField
    (
        const dgGeomMesh* dgMesh,
        const cellDof<Type>* dof
    );

    //- Copy constructor
    cellGaussField(const cellGaussField<Type>& other);

    //- Construct with uniform initial value
    cellGaussField
    (
        const label cellID,
        const dgGeomMesh* dgMesh,
        const Type& initVal
    );

    //- Construct with zero-initialised values
    cellGaussField
    (
        const label cellID,
        const dgGeomMesh* dgMesh
    );

    // Access -----------------------------------------------------------------

    inline Type& operator[](const label gp)
    {
        return values_[gp];
    }

    inline const Type& operator[](const label gp) const
    {
        return values_[gp];
    }

    inline label size() const
    {
        return nGauss_;
    }

    inline label cellID() const
    {
        return cellID_;
    }

    inline const dgGeomMesh* dgMesh() const
    {
        if (!dgMesh_)
        {
            FatalErrorInFunction << "dgMesh_ is null" << abort(FatalError);
        }
        return dgMesh_;
    }

    //- Interpolate DOF → Gauss values
    void interpolateFromDof();

    // Operators --------------------------------------------------------------

    //- Copy assignment
    cellGaussField<Type>& operator=
    (
        const cellGaussField<Type>& other
    );

    //- Assign constant value to all Gauss points
    cellGaussField<Type>& operator=
    (
        const Type& value
    );

    // Printing ---------------------------------------------------------------

    friend Ostream& operator<<
    (
        Ostream& os,
        const cellGaussField<Type>& cgf
    )
    {
        os << "cellGaussField for cell " << cgf.cellID() << ":\n";
        for (label i = 0; i < cgf.size(); ++i)
        {
            os << "  Gauss point " << i << ": " << cgf[i] << nl;
        }
        return os;
    }

    // Scalar/vector specialisations declared here
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

    // Limiting
    void limitMax(const scalar maxVal);
    void limitMin(const scalar minVal);
    void limitRange(const scalar minVal, const scalar maxVal);
};

} // End namespace Foam

// Template implementation
#include "cellGaussFieldI.H"

// Specialisations
#include "cellGaussFieldScalarI.H"
#include "cellGaussFieldVectorI.H"

#endif

// ************************************************************************* //

