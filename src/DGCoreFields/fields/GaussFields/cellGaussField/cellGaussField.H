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
    Foam::cellGaussField

Description
    Gauss-point field values for a DG cell.

    Stores field values evaluated at Gauss integration points
    within a single DG cell.

SourceFiles
    cellGaussFieldI.H
    cellGaussFieldScalarI.H
    cellGaussFieldVectorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellGaussField_H
#define Foam_cellGaussField_H

#include "label.H"
#include "List.H"
#include "refCount.H"
#include "dgGeomMesh.H"
#include "dgGeomCell.H"
#include "cellDof.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cellGaussField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::cellGaussField
 * \brief Field values at Gauss points of a DG cell
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * cellGaussField stores field values evaluated at Gauss integration
 * points inside a single DG cell.
 *
 * The values are typically obtained by interpolating cell-local
 * degrees of freedom (DoF) using basis functions provided by the
 * DG geometric cell.
 *
 * This class does not own geometric data or DoF storage; it holds
 * non-owning pointers to dgGeomMesh, dgGeomCell, and cellDof<Type>.
 */
template<class Type>
class cellGaussField
:
    public refCount
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** \brief Mesh cell index associated with this Gauss field. */
        label cellID_;

        /** \brief Non-owning pointer to the DG geometric mesh. */
        const dgGeomMesh* dgMesh_;

        /** \brief Non-owning pointer to the DG geometric cell. */
        const dgGeomCell* cell_;

        /** \brief Number of interior Gauss points in the cell. */
        label nGauss_;

        /** \brief Number of local degrees of freedom in the cell. */
        label nDof_;

        /** \brief Non-owning pointer to the cell-local DoF container. */
        const cellDof<Type>* dof_;

        /** \brief Field values stored at the cell Gauss integration points. */
        List<Type> values_;


public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Default constructor
         *
         * Required for List containers.
         */
        cellGaussField();

        /**
         * \brief Construct from DG mesh and DoF pointer
         *
         * Initializes Gauss-point storage based on the cell geometry
         * and DoF information.
         *
         * \param dgMesh DG geometric mesh containing the cell.
         * \param dof Cell-local DoF data used to size and later interpolate the field.
         */
        cellGaussField
        (
            const dgGeomMesh* dgMesh,
            const cellDof<Type>* dof
        );

        /**
         * \brief Copy constructor.
         *
         * \param other Source cell Gauss field to copy from.
         */
        cellGaussField(const cellGaussField<Type>& other);

        /**
         * \brief Construct with a uniform initial value.
         *
         * \param cellID Mesh cell index.
         * \param dgMesh DG geometric mesh containing the cell.
         * \param initVal Initial value assigned at every Gauss point.
         */
        cellGaussField
        (
            const label cellID,
            const dgGeomMesh* dgMesh,
            const Type& initVal
        );

        /**
         * \brief Construct with zero-initialized values.
         *
         * \param cellID Mesh cell index.
         * \param dgMesh DG geometric mesh containing the cell.
         */
        cellGaussField
        (
            const label cellID,
            const dgGeomMesh* dgMesh
        );

        /**
         * \brief Move constructor
         *
         * Transfers internal storage without copying.
         *
         * \param other Source cell Gauss field whose storage is moved from.
         */
        cellGaussField(cellGaussField<Type>&& other) noexcept;

        /**
         * \brief Construct by consuming a temporary cellGaussField.
         *
         * \param t Temporary cell Gauss field whose storage is consumed.
         */
        cellGaussField
        (
            const tmp<cellGaussField<Type>>& t
        );

    //====================================================================//
    //                              Access
    //====================================================================//

        /**
         * \brief Return writable access to one cell Gauss-point value.
         *
         * \param gp Gauss-point index inside the cell.
         * \return Reference to the stored value.
         */
        inline Type& operator[](const label gp)
        {
            return values_[gp];
        }

        /**
         * \brief Return read-only access to one cell Gauss-point value.
         *
         * \param gp Gauss-point index inside the cell.
         * \return Const reference to the stored value.
         */
        inline const Type& operator[](const label gp) const
        {
            return values_[gp];
        }

        /** \brief Return the number of Gauss points stored for the cell. */
        inline label size() const
        {
            return nGauss_;
        }

        /** \brief Return the associated mesh cell index. */
        inline label cellID() const
        {
            return cellID_;
        }

        /**
         * \brief Return DG geometric mesh
         *
         * Fatal error if pointer is null.
         *
         * \return Pointer to the DG geometric mesh.
         */
        inline const dgGeomMesh* dgMesh() const
        {
            if (!dgMesh_)
            {
                FatalErrorInFunction << "dgMesh_ is null" << abort(FatalError);
            }
            return dgMesh_;
        }

        /**
         * \brief Return DG geometric cell
         *
         * Fatal error if pointer is null.
         *
         * \return Pointer to the DG geometric cell.
         */
        inline const dgGeomCell* cell() const
        {
            if (!cell_)
            {
                FatalErrorInFunction << "cell_ is null" << abort(FatalError);
            }
            return cell_;
        }

    //====================================================================//
    //                        Interpolation
    //====================================================================//

        /**
         * \brief Interpolate DoF values to Gauss points
         *
         * Evaluates field values at Gauss points using basis functions
         * and the associated DoF values.
         */
        void interpolateFromDof();

    //====================================================================//
    //                           Operators
    //====================================================================//

        /**
         * \brief Copy assignment operator.
         *
         * \param other Source cell Gauss field to copy from.
         * \return Reference to this object.
         */
        cellGaussField<Type>& operator=
        (
            const cellGaussField<Type>& other
        );

        /**
         * \brief Assign a uniform value to all Gauss points.
         *
         * \param value Uniform value assigned at every Gauss point.
         * \return Reference to this object.
         */
        cellGaussField<Type>& operator=
        (
            const Type& value
        );

        /**
         * \brief Move assignment
         *
         * Transfers internal storage without copying.
         *
         * \param other Source cell Gauss field whose storage is moved from.
         * \return Reference to this object.
         */
        cellGaussField<Type>& operator=
        (
            cellGaussField<Type>&& other
        ) noexcept;

        /**
         * \brief Assign by consuming a temporary cellGaussField.
         *
         * \param t Temporary cell Gauss field whose storage is consumed.
         * \return Reference to this object.
         */
        cellGaussField<Type>& operator=
        (
            const tmp<cellGaussField<Type>>& t
        );

    //====================================================================//
    //                             Output
    //====================================================================//

        /**
         * \brief Write the cell Gauss-point values to an output stream.
         *
         * \param os Output stream.
         * \param cgf Cell Gauss field to serialize.
         * \return Reference to the output stream.
         */
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

    //====================================================================//
    //                    Component access (vector fields)
    //====================================================================//

        /** \brief Return the x-component as a scalar cellGaussField. */
        cellGaussField<scalar> x() const;

        /** \brief Return the y-component as a scalar cellGaussField. */
        cellGaussField<scalar> y() const;

        /** \brief Return the z-component as a scalar cellGaussField. */
        cellGaussField<scalar> z() const;

        /** \brief Replace the x-component using a scalar cellGaussField. */
        void setX(const cellGaussField<scalar>& fx);

        /** \brief Replace the y-component using a scalar cellGaussField. */
        void setY(const cellGaussField<scalar>& fy);

        /** \brief Replace the z-component using a scalar cellGaussField. */
        void setZ(const cellGaussField<scalar>& fz);

        /** \brief Replace all vector components using scalar cellGaussFields. */
        void setXYZ
        (
            const cellGaussField<scalar>& fx,
            const cellGaussField<scalar>& fy,
            const cellGaussField<scalar>& fz
        );

    //====================================================================//
    //                              Limiting
    //====================================================================//

        /** \brief Limit all stored values to a maximum value. */
        void limitMax(const scalar maxVal);

        /** \brief Limit all stored values to a minimum value. */
        void limitMin(const scalar minVal);

        /** \brief Limit all stored values to the range [minVal, maxVal]. */
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
