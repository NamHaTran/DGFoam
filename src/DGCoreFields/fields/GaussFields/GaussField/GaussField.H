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
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGFoam is distributed in the hope that it will be useful, but WITHOUT
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
    GaussFieldI.H
    GaussFieldVectorI.H
    GaussFieldScalarI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_GaussField_H
#define Foam_GaussField_H

#include "cellGaussField.H"
#include "faceGaussField.H"
#include "cellDof.H"
#include "dofField.H"
#include "dgGeomMesh.H"
#include "refCount.H"

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
:
    public refCount
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** \brief Non-owning pointer to the source DoF field, when available. */
        const dofField<Type>* dofField_;

        /** \brief Mesh cell index associated with this Gauss field. */
        label cellID_;

        /** \brief Non-owning pointer to the DG geometric mesh. */
        const dgGeomMesh* mesh_;

        /** \brief Cell-interior Gauss-point values. */
        cellGaussField<Type> cellField_;

        /** \brief Face Gauss-point values on the cell boundary. */
        faceGaussField<Type> faceField_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** \brief Construct an empty GaussField. */
        GaussField();

        /**
         * \brief Construct from DoF field
         *
         * The Gauss values are initialized by interpolating
         * the DOF values of the given cell.
         *
         * \param dofField Source DoF field.
         * \param cellID Mesh cell index.
         * \param mesh DG geometric mesh containing the cell.
         */
        GaussField
        (
            const dofField<Type>* dofField,
            label cellID,
            const dgGeomMesh* mesh
        );

        /**
         * \brief Construct an empty GaussField for one cell.
         *
         * \param cellID Mesh cell index.
         * \param mesh DG geometric mesh containing the cell.
         */
        GaussField
        (
            label cellID,
            const dgGeomMesh* mesh
        );

        /**
         * \brief Construct a GaussField with a uniform initial value.
         *
         * \param cellID Mesh cell index.
         * \param mesh DG geometric mesh containing the cell.
         * \param initVal Initial value assigned to all cell and face Gauss points.
         */
        GaussField
        (
            label cellID,
            const dgGeomMesh* mesh,
            const Type& initVal
        );

        /**
         * \brief Copy constructor.
         *
         * \param other Source GaussField to copy from.
         */
        GaussField(const GaussField<Type>& other);

        /**
         * \brief Construct by consuming a temporary GaussField.
         *
         * \param t Temporary GaussField whose storage is consumed.
         */
        GaussField
        (
            const tmp<GaussField<Type>>& t
        );

        /**
         * \brief Construct a temporary GaussField for one cell.
         *
         * \param cellID Mesh cell index.
         * \param mesh DG geometric mesh containing the cell.
         * \return Temporary GaussField.
         */
        static tmp<GaussField<Type>> New
        (
            label cellID,
            const dgGeomMesh* mesh
        );

        /**
         * \brief Move constructor.
         *
         * \param other Source GaussField whose storage is moved from.
         */
        GaussField
        (
            GaussField<Type>&& other
        ) noexcept;


    //====================================================================//
    //                           Assignment
    //====================================================================//

        /**
         * \brief Copy assignment operator.
         *
         * \param other Source GaussField to copy from.
         * \return Reference to this object.
         */
        GaussField<Type>& operator=(const GaussField<Type>& other);

        /**
         * \brief Assign a uniform value to all stored Gauss points.
         *
         * \param value Uniform value assigned to cell and face Gauss points.
         * \return Reference to this object.
         */
        GaussField<Type>& operator=(const Type& value);

        /**
         * \brief Move assignment operator.
         *
         * \param other Source GaussField whose storage is moved from.
         * \return Reference to this object.
         */
        GaussField<Type>& operator=
        (
            GaussField<Type>&& other
        ) noexcept;

        /**
         * \brief Assign by consuming a temporary GaussField.
         *
         * \param t Temporary GaussField whose storage is consumed.
         * \return Reference to this object.
         */
        GaussField<Type>& operator=
        (
            const tmp<GaussField<Type>>& t
        );

    //====================================================================//
    //                             Access
    //====================================================================//

        /**
         * \brief Return DG geometric mesh
         *
         * Fatal error if mesh pointer is null.
         *
         * \return Pointer to the DG geometric mesh.
         */
        inline const dgGeomMesh* dgMesh() const;

        /** \brief Return the associated mesh cell index. */
        inline label cellID() const { return cellID_; }

        /** \brief Return the pointer to the source DoF field, if any. */
        inline const dofField<Type>* dofFieldPtr() const { return dofField_; }

        /** \brief Return writable access to the cell-interior Gauss field. */
        inline cellGaussField<Type>& cellField()       { return cellField_; }

        /** \brief Return read-only access to the cell-interior Gauss field. */
        inline const cellGaussField<Type>& cellField() const { return cellField_; }

        /** \brief Return writable access to the face Gauss field. */
        inline faceGaussField<Type>& faceField()       { return faceField_; }

        /** \brief Return read-only access to the face Gauss field. */
        inline const faceGaussField<Type>& faceField() const { return faceField_; }

        /**
         * \brief Access neighbour DoFs for one local face of this cell.
         *
         * Internal faces return directly linked neighbour-cell DoFs.
         * Processor faces return the most recently synchronized buffer.
         *
         * \param localFaceI Local face index.
         * \return Read-only neighbour cell DoF data.
         */
        inline const cellDof<Type>& neighbourCellDof
        (
            const label localFaceI
        ) const
        {
            return faceField_.neighbourCellDof(localFaceI);
        }


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

        /** \brief Return writable access to one cell-interior Gauss value. */
        inline Type& cellValueAt(label gpI);

        /** \brief Return read-only access to one cell-interior Gauss value. */
        inline const Type& cellValueAt(label gpI) const;

        /** \brief Return writable access to one plus-side face Gauss value. */
        inline Type& facePlusValueAt(label faceI, label gpI);

        /** \brief Return writable access to one minus-side face Gauss value. */
        inline Type& faceMinusValueAt(label faceI, label gpI);

        /** \brief Return read-only access to one plus-side face Gauss value. */
        inline const Type& facePlusValueAt(label faceI, label gpI) const;

        /** \brief Return read-only access to one minus-side face Gauss value. */
        inline const Type& faceMinusValueAt(label faceI, label gpI) const;


    //====================================================================//
    //                             Output
    //====================================================================//

        /**
         * \brief Write the combined cell and face Gauss values to a stream.
         *
         * \param os Output stream.
         * \param gf Gauss field to serialize.
         * \return Reference to the output stream.
         */
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

        /** \brief Return the x-component as a scalar GaussField. */
        GaussField<scalar> x() const;

        /** \brief Return the y-component as a scalar GaussField. */
        GaussField<scalar> y() const;

        /** \brief Return the z-component as a scalar GaussField. */
        GaussField<scalar> z() const;

        /** \brief Replace the x-component using a scalar GaussField. */
        void setX(const GaussField<scalar>& fx);

        /** \brief Replace the y-component using a scalar GaussField. */
        void setY(const GaussField<scalar>& fy);

        /** \brief Replace the z-component using a scalar GaussField. */
        void setZ(const GaussField<scalar>& fz);

        /** \brief Replace all vector components using scalar GaussFields. */
        void setXYZ
        (
            const GaussField<scalar>& fx,
            const GaussField<scalar>& fy,
            const GaussField<scalar>& fz
        );

    //====================================================================//
    //                              Limiting
    //====================================================================//

        /** \brief Limit all stored values to a maximum value. */
        void limitMax(scalar maxVal);

        /** \brief Limit all stored values to a minimum value. */
        void limitMin(scalar minVal);

        /** \brief Limit all stored values to the range [minVal, maxVal]. */
        void limitRange(scalar minVal, scalar maxVal);
};

} // End namespace Foam


// * * * * * * * * Template implementation * * * * * * * * //
#include "GaussFieldI.H"
#include "GaussFieldVectorI.H"
#include "GaussFieldScalarI.H"

#endif
// ************************************************************************* //
