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
    Foam::dofField<Type>

Description
    Degree-of-freedom field for DG discretization.

    Stores all cell-wise degrees of freedom and provides a bridge
    between DG DoF storage and OpenFOAM volField representation.

SourceFiles
    dofFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dofField_H
#define Foam_dofField_H

#include "volFields.H"
#include "cellDof.H"
#include "dgGeomMesh.H"
#include "IOobject.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dofField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dofField
 * \brief DG degree-of-freedom field
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * dofField stores degrees of freedom (DoFs) for all cells in a DG mesh.
 *
 * The class maintains:
 * - Cell-local DoF storage via cellDof<Type>
 * - A set of OpenFOAM volFields, one per DoF index
 *
 * The OpenFOAM volFields allow DG DoF data to be written to disk,
 * visualized, or post-processed using standard OpenFOAM tools.
 *
 * The class does not perform any numerical operations; it only
 * manages data storage and synchronization between DG and FV views.
 */
template<class Type>
class dofField
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** \brief DG geometric mesh that defines the cell-wise DoF layout. */
        const dgGeomMesh& dgMesh_;

        /** \brief Cell-local degree-of-freedom storage for the full mesh. */
        List<cellDof<Type>> cellDofs_;

        /** \brief Number of cells represented in the field. */
        label nCells_;

        /**
         * OpenFOAM volFields corresponding to each DoF index.
         *
         * Each volField stores one DoF component over all cells.
         */
        List<autoPtr<GeometricField<Type, fvPatchField, volMesh>>> foamFields_;

        /** \brief Number of local degrees of freedom for each cell. */
        List<label> dofPerCell_;

        /** \brief Maximum DoF count over all cells in the mesh. */
        label maxDoF_;

        /**
         * \brief Refresh boundary faces from the adjacent cell values.
         *
         * This keeps the on-disk boundaryField section synchronized with the
         * current cell values before a volField is written.
         *
         * \param dofI DoF index whose boundary values should be refreshed.
         */
        void updateBoundaryField(const label dofI);

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct a DoF field from OpenFOAM I/O information.
         *
         * Creates one OpenFOAM volume field per DoF index, reads available
         * field files from disk, and reconstructs the cell-local DG DoF data.
         *
         * \param name Base field name.
         * \param instance Time or instance directory.
         * \param rOpt Read option used for the field object.
         * \param wOpt Write option used for the generated volFields.
         * \param mesh DG geometric mesh associated with the field.
         */
        dofField
        (
            const word& name,
            const fileName& instance,
            IOobject::readOption rOpt,
            IOobject::writeOption wOpt,
            const dgGeomMesh& mesh
        );

        /**
         * \brief Construct a zero-initialized DoF field.
         *
         * Allocates cell-local DoF storage for the given DG mesh and
         * initializes all degrees of freedom to zero.
         *
         * \param mesh DG geometric mesh associated with the field.
         */
        dofField
        (
            const dgGeomMesh& mesh
        );

    //====================================================================//
    //                              Access
    //====================================================================//

        /** \brief Return read-only access to the full list of cell DoFs. */
        inline const List<cellDof<Type>>& cellDofList() const
        {
            return cellDofs_;
        }

        /** \brief Return writable access to the full list of cell DoFs. */
        inline List<cellDof<Type>>& cellDofList()
        {
            return cellDofs_;
        }

        /**
         * \brief Return writable access to the DoFs of one cell.
         *
         * \param i Cell index.
         * \return Reference to the cell-local DoF container.
         */
        inline cellDof<Type>& operator[](const label i)
        {
            return cellDofs_[i];
        }

        /**
         * \brief Return read-only access to the DoFs of one cell.
         *
         * \param i Cell index.
         * \return Const reference to the cell-local DoF container.
         */
        inline const cellDof<Type>& operator[](const label i) const
        {
            return cellDofs_[i];
        }

    //====================================================================//
    //                         Synchronization
    //====================================================================//

        /**
         * \brief Update the OpenFOAM volFields from one cell's DG DoFs.
         *
         * Copies the DoFs of a single DG cell into the corresponding
         * OpenFOAM volFields. This is the preferred hot-path synchronization
         * helper during explicit stage updates.
         *
         * \param cellID Cell index whose DoFs are copied to the volFields.
         */
        void updateCellDof(const label cellID);

        /**
         * \brief Write all OpenFOAM DoF volFields to disk.
         */
        void write() const;
};

} // End namespace Foam

#include "dofFieldI.H"

#endif

// ************************************************************************* //
