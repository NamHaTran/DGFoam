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
    dofField.C

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

        /** Reference to DG geometric mesh */
        const dgGeomMesh& dgMesh_;

        /** Cell-local DoF storage */
        List<cellDof<Type>> cellDofs_;

        /** Number of cells in the mesh */
        label nCells_;

        /**
         * OpenFOAM volFields corresponding to each DoF index.
         *
         * Each volField stores one DoF component over all cells.
         */
        List<autoPtr<GeometricField<Type, fvPatchField, volMesh>>> foamFields_;

        /** Number of DoFs per cell */
        List<label> dofPerCell_;

        /** Maximum number of DoFs across all cells */
        label maxDoF_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct from field name and I/O options
         *
         * Creates OpenFOAM volFields corresponding to DG DoFs
         * and reads data according to the provided IOobject options.
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
         * \brief Construct zero-initialized DoF field
         *
         * Creates an empty DoF field associated with the given
         * DG geometric mesh.
         */
        dofField
        (
            const dgGeomMesh& mesh
        );

    //====================================================================//
    //                              Access
    //====================================================================//

        /** Return const access to cell DoF list */
        inline const List<cellDof<Type>>& cellDofList() const
        {
            return cellDofs_;
        }

        /** Return non-const access to cell DoF list */
        inline List<cellDof<Type>>& cellDofList()
        {
            return cellDofs_;
        }

        /** Return non-const access to DoFs of a given cell */
        inline cellDof<Type>& operator[](const label i)
        {
            return cellDofs_[i];
        }

        /** Return const access to DoFs of a given cell */
        inline const cellDof<Type>& operator[](const label i) const
        {
            return cellDofs_[i];
        }

    //====================================================================//
    //                         Synchronization
    //====================================================================//

        /**
         * \brief Update OpenFOAM volFields from DG DoF values
         *
         * Copies cell-local DG DoF values into the corresponding
         * OpenFOAM volFields.
         */
        void updateFoamFields();

        /**
         * \brief Write OpenFOAM volFields to disk
         */
        void write() const;
};

} // End namespace Foam

#include "dofFieldI.H"

#endif

// ************************************************************************* //
