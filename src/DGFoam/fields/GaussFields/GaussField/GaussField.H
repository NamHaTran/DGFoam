/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Foam::GaussField

Description
    Wrapper class for Gauss-point field data of a DG cell.

    This class aggregates Gauss-point values inside a cell and on its faces.
    It provides a unified interface to access and operate on:
      - Cell interior Gauss-point values
      - Face Gauss-point values (minus/plus sides)

SourceFiles
    GaussField.C

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

/**
 * \class Foam::GaussField
 * \brief Combined Gauss-point field for a DG cell
 *
 * \tparam Type Field value type (scalar, vector, tensor, etc.)
 *
 * \details
 * GaussField acts as a high-level container that combines:
 *  - cellGaussField<Type> for cell-interior Gauss points
 *  - faceGaussField<Type> for Gauss points on cell faces
 *
 * The class can be constructed:
 *  - From a dofField<Type> to interpolate DOF values to Gauss points
 *  - As an empty container to be filled later
 *  - With a uniform initial value
 *
 * GaussField does not own mesh or DoF storage. All pointers are non-owning.
 * It is intended to be stored per-cell inside higher-level DG field classes.
 */
template<class Type>
class GaussField
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Pointer to DoF field (read-only, non-owning) */
        const dofField<Type>* dofField_;

        /** Cell index */
        label cellID_;

        /** Pointer to DG geometric mesh (non-owning) */
        const dgGeomMesh* mesh_;

        /** Gauss-point field inside the cell */
        cellGaussField<Type> cellField_;

        /** Gauss-point field on cell faces */
        faceGaussField<Type> faceField_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** Default constructor */
        GaussField();

        /**
         * \brief Construct from DoF field
         *
         * The Gauss values are initialized by interpolating
         * the DOF values of the given cell.
         */
        GaussField
        (
            const dofField<Type>* dofField,
            label cellID,
            const dgGeomMesh* mesh
        );

        /**
         * \brief Construct empty GaussField for a cell
         */
        GaussField
        (
            label cellID,
            const dgGeomMesh* mesh
        );

        /**
         * \brief Construct GaussField with uniform initial value
         */
        GaussField
        (
            label cellID,
            const dgGeomMesh* mesh,
            const Type& initVal
        );

        /** Copy constructor */
        GaussField(const GaussField<Type>& other);


    //====================================================================//
    //                           Assignment
    //====================================================================//

        /** Copy assignment */
        GaussField<Type>& operator=(const GaussField<Type>& other);

        /** Assign uniform value */
        GaussField<Type>& operator=(const Type& value);


    //====================================================================//
    //                             Access
    //====================================================================//

        /**
         * \brief Return DG geometric mesh
         *
         * Fatal error if mesh pointer is null.
         */
        inline const dgGeomMesh* dgMesh() const;

        /** Return cell index */
        inline label cellID() const { return cellID_; }

        /** Return pointer to DoF field */
        inline const dofField<Type>* dofFieldPtr() const { return dofField_; }

        /** Access cell Gauss field */
        inline cellGaussField<Type>& cellField()       { return cellField_; }

        /** Const access to cell Gauss field */
        inline const cellGaussField<Type>& cellField() const { return cellField_; }

        /** Access face Gauss field */
        inline faceGaussField<Type>& faceField()       { return faceField_; }

        /** Const access to face Gauss field */
        inline const faceGaussField<Type>& faceField() const { return faceField_; }


    //====================================================================//
    //                       High-level operations
    //====================================================================//

        /**
         * \brief Interpolate DOF values to Gauss points
         *
         * Updates both cell and face Gauss fields using
         * the associated dofField<Type>.
         */
        void interpolateFromDof();


    //====================================================================//
    //                      Direct Gauss access
    //====================================================================//

        /** Access cell Gauss value by index */
        inline Type& cellValueAt(label gpI);

        /** Const access to cell Gauss value */
        inline const Type& cellValueAt(label gpI) const;

        /** Access plus-side face Gauss value */
        inline Type& facePlusValueAt(label faceI, label gpI);

        /** Access minus-side face Gauss value */
        inline Type& faceMinusValueAt(label faceI, label gpI);

        /** Const access to plus-side face Gauss value */
        inline const Type& facePlusValueAt(label faceI, label gpI) const;

        /** Const access to minus-side face Gauss value */
        inline const Type& faceMinusValueAt(label faceI, label gpI) const;


    //====================================================================//
    //                             Output
    //====================================================================//

        /** Stream output */
        friend Ostream& operator<<(Ostream& os, const GaussField& gf)
        {
            os << "GaussField for cell " << gf.cellID() << ":\n";
            os << "Cell Gauss Field:\n" << gf.cellField() << "\n";
            os << "Face Gauss Field:\n" << gf.faceField() << "\n";
            return os;
        }


    //====================================================================//
    //                Component access (vector fields)
    //====================================================================//

        /** Return x-component */
        GaussField<scalar> x() const;

        /** Return y-component */
        GaussField<scalar> y() const;

        /** Return z-component */
        GaussField<scalar> z() const;

        /** Set x-component */
        void setX(const GaussField<scalar>& fx);

        /** Set y-component */
        void setY(const GaussField<scalar>& fy);

        /** Set z-component */
        void setZ(const GaussField<scalar>& fz);

        /** Set all components */
        void setXYZ
        (
            const GaussField<scalar>& fx,
            const GaussField<scalar>& fy,
            const GaussField<scalar>& fz
        );

    //====================================================================//
    //                              Limiting
    //====================================================================//

        /** Limit maximum value */
        void limitMax(scalar maxVal);

        /** Limit minimum value */
        void limitMin(scalar minVal);

        /** Limit value range */
        void limitRange(scalar minVal, scalar maxVal);
};

} // End namespace Foam


// * * * * * * * * Template implementation * * * * * * * * //
#include "GaussFieldI.H"
#include "GaussFieldVectorI.H"
#include "GaussFieldScalarI.H"

#endif
// ************************************************************************* //
