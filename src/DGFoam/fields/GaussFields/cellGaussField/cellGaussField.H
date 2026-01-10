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
    Foam::cellGaussField

Description
    Gauss-point field values for a DG cell.

    Stores field values evaluated at Gauss integration points
    within a single DG cell.

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
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Cell index in the mesh */
        label cellID_;

        /** Pointer to DG geometric mesh (non-owning) */
        const dgGeomMesh* dgMesh_;

        /** Pointer to DG geometric cell (non-owning) */
        const dgGeomCell* cell_;

        /** Number of Gauss points in this cell */
        label nGauss_;

        /** Number of degrees of freedom in this cell */
        label nDof_;

        /** Pointer to cell-local DoF container (non-owning) */
        const cellDof<Type>* dof_;

        /** Values at Gauss integration points */
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
         */
        cellGaussField
        (
            const dgGeomMesh* dgMesh,
            const cellDof<Type>* dof
        );

        /**
         * \brief Copy constructor
         */
        cellGaussField(const cellGaussField<Type>& other);

        /**
         * \brief Construct with uniform initial value
         */
        cellGaussField
        (
            const label cellID,
            const dgGeomMesh* dgMesh,
            const Type& initVal
        );

        /**
         * \brief Construct with zero-initialized values
         */
        cellGaussField
        (
            const label cellID,
            const dgGeomMesh* dgMesh
        );

    //====================================================================//
    //                              Access
    //====================================================================//

        /** Non-const access to Gauss-point value */
        inline Type& operator[](const label gp)
        {
            return values_[gp];
        }

        /** Const access to Gauss-point value */
        inline const Type& operator[](const label gp) const
        {
            return values_[gp];
        }

        /** Return number of Gauss points */
        inline label size() const
        {
            return nGauss_;
        }

        /** Return cell index */
        inline label cellID() const
        {
            return cellID_;
        }

        /**
         * \brief Return DG geometric mesh
         *
         * Fatal error if pointer is null.
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
         * \brief Copy assignment
         */
        cellGaussField<Type>& operator=
        (
            const cellGaussField<Type>& other
        );

        /**
         * \brief Assign uniform value to all Gauss points
         */
        cellGaussField<Type>& operator=
        (
            const Type& value
        );

    //====================================================================//
    //                             Output
    //====================================================================//

        /**
         * \brief Stream output of Gauss-point values
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

        /** Return x-component as scalar Gauss field */
        cellGaussField<scalar> x() const;

        /** Return y-component as scalar Gauss field */
        cellGaussField<scalar> y() const;

        /** Return z-component as scalar Gauss field */
        cellGaussField<scalar> z() const;

        /** Set x-component from scalar Gauss field */
        void setX(const cellGaussField<scalar>& fx);

        /** Set y-component from scalar Gauss field */
        void setY(const cellGaussField<scalar>& fy);

        /** Set z-component from scalar Gauss field */
        void setZ(const cellGaussField<scalar>& fz);

        /** Set all components from scalar Gauss fields */
        void setXYZ
        (
            const cellGaussField<scalar>& fx,
            const cellGaussField<scalar>& fy,
            const cellGaussField<scalar>& fz
        );

    //====================================================================//
    //                              Limiting
    //====================================================================//

        /** Limit maximum value */
        void limitMax(const scalar maxVal);

        /** Limit minimum value */
        void limitMin(const scalar minVal);

        /** Limit value range */
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
