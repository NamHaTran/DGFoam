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
    Foam::dgField

Description
    Discontinuous Galerkin field container.

    Stores degrees of freedom (DoF) and corresponding GaussField
    representations for all cells in the mesh.

SourceFiles
    dgField.C
    dgFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgField_H
#define Foam_dgField_H

#include "regIOobject.H"
#include "cellDof.H"
#include "dofField.H"
#include "GaussField.H"
#include "dgGeomMesh.H"
#include "dgGeneralBoundaryManager.H"
#include "tensor.H"

namespace Foam
{

template<class Type>
class dgProcessorBoundaryManager;

/*---------------------------------------------------------------------------*\
                          Class dgField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgField
 * \brief DG field storing DoFs and Gauss-point values
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * dgField represents a field discretized using the Discontinuous
 * Galerkin method over the entire mesh.
 *
 * The class stores:
 * - Degrees of freedom (DoF) per cell via dofField<Type>
 * - Precomputed GaussField<Type> values for each cell
 *
 * A dgField may own a DoF field (primary DG variables) or act
 * as an auxiliary field without DoFs (derived quantities).
 *
 * The class derives from regIOobject to integrate with OpenFOAM
 * I/O and objectRegistry mechanisms.
 */
template<class Type>
class dgField
:
    public regIOobject
{
    private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** \brief DG geometric mesh on which the field is defined. */
        const dgGeomMesh& mesh_;

        /** \brief Flag indicating whether the field is constructed with DoF storage. */
        bool hasDof_;

        /** \brief Optional DoF storage, null when the field is Gauss-value only. */
        autoPtr<dofField<Type>> dofPtr_;

        /** \brief Per-cell Gauss-point representation of the DG field. */
        List<GaussField<Type>> gaussFields_;

        /**
         * Lazy processor-boundary synchronizer for this field.
         *
         * Allocated on first synch() call so call sites can simply use
         * field.synch() without wiring an external manager everywhere.
         */
        autoPtr<dgProcessorBoundaryManager<Type>> procSyncPtr_;

        /**
         * \brief Build a temporary volField containing cell-average values.
         *
         * The returned field is used purely as an I/O adapter so DG fields can
         * be written in a ParaView-readable OpenFOAM volume-field format.
         */
        tmp<GeometricField<Type, fvPatchField, volMesh>> cellAverageVolField() const;

public:

    //- Runtime type information
    TypeName("dgField");

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct a DG field on a DG geometric mesh.
         *
         * Optionally allocates DoF storage and initializes the corresponding
         * Gauss-field representation for every cell in the mesh.
         *
         * \param io OpenFOAM I/O descriptor for the field object.
         * \param mesh DG geometric mesh on which the field lives.
         * \param hasDof Set to \c true to allocate a DoF field, or \c false
         *     for Gauss-value-only fields.
         */
        dgField
        (
            const IOobject& io,
            const dgGeomMesh& mesh,
            bool hasDof = true
        );

        /** \brief Copy construction is disabled. */
        dgField(const dgField&) = delete;

        /** \brief Copy assignment is disabled. */
        dgField& operator=(const dgField&) = delete;

    //====================================================================//
    //                           Access
    //====================================================================//

        /** \brief Return true if this field currently owns DoF storage. */
        inline bool hasDof() const
        {
            return dofPtr_ != nullptr;
        }

        /**
         * \brief Return writable access to the DoF field.
         *
         * Fatal error if this field does not own DoF data.
         *
         * \return Reference to the owned DoF field.
         */
        inline dofField<Type>& dof()
        {
            if (!dofPtr_)
            {
                FatalErrorInFunction
                    << "Attempt to access non-existing DOF field"
                    << abort(FatalError);
            }
            return *dofPtr_;
        }

        /**
         * \brief Return read-only access to the DoF field.
         *
         * Fatal error if this field does not own DoF data.
         *
         * \return Const reference to the owned DoF field.
         */
        inline const dofField<Type>& dof() const
        {
            if (!dofPtr_)
            {
                FatalErrorInFunction
                    << "Attempt to access non-existing DOF field"
                    << abort(FatalError);
            }
            return *dofPtr_;
        }

        /** \brief Return writable access to all per-cell Gauss fields. */
        inline List<GaussField<Type>>& gaussFields()
        {
            return gaussFields_;
        }

        /** \brief Return read-only access to all per-cell Gauss fields. */
        inline const List<GaussField<Type>>& gaussFields() const
        {
            return gaussFields_;
        }

        /** \brief Return the number of cells covered by the DG field. */
        inline label nCells() const
        {
            return gaussFields_.size();
        }

        /** \brief Return the registered object name of the field. */
        inline const word& name() const
        {
            return regIOobject::name();
        }

        /**
         * \brief Synchronize processor-boundary data for this DG field.
         *
         * On parallel runs this updates processor-face plus states and
         * buffered neighbour DoFs so downstream DG operators can consume a
         * consistent interface state via a simple field-local API.
         */
        void synch();

        /**
         * \brief Refresh all Gauss values from the owned DoF field.
         *
         * Cell-interior values are interpolated for every cell.  For local
         * internal faces, only the owner cell reconstructs the face traces;
         * the neighbour cell copies the owner plus/minus traces with sides
         * swapped.  Boundary and processor faces keep the existing local
         * interpolation behaviour.
         */
        void interpolateFromDof();

        /**
         * \brief Access neighbour DoFs for one face of one cell.
         *
         * This helper hides whether the neighbour is local (internal face)
         * or comes from a processor-face synchronization buffer.
         *
         * \param cellID Owner cell index.
         * \param localFaceI Local face index inside the owner cell.
         * \return Read-only neighbour cell DoF data for that face.
         */
        inline const cellDof<Type>& neighbourCellDof
        (
            const label cellID,
            const label localFaceI
        ) const
        {
            if (cellID < 0 || cellID >= gaussFields_.size())
            {
                FatalErrorInFunction
                    << "Cell index " << cellID << " is out of range [0, "
                    << gaussFields_.size() - 1 << ']'
                    << abort(FatalError);
            }

            return gaussFields_[cellID].neighbourCellDof(localFaceI);
        }

    //====================================================================//
    //                              I/O
    //====================================================================//

        /**
         * \brief Write the field data payload to a stream.
         *
         * \param os Output stream receiving the serialized data.
         * \return True if the write succeeds.
         */
        bool writeData(Ostream& os) const;

        /**
         * \brief Write this DG field as a standard OpenFOAM volField.
         *
         * This overrides regIOobject writing so the on-disk file header matches
         * `volScalarField`/`volVectorField` etc, which ParaView can read.
         */
        bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const override;

    //====================================================================//
    //                           Assignment
    //====================================================================//

        /**
         * \brief Assign a uniform value to all Gauss fields.
         *
         * \param value Uniform value assigned at every stored Gauss location.
         * \return Reference to this field.
         */
        dgField<Type>& operator=(const Type& value);
};

} // End namespace Foam

#include "dgFieldI.H"

#endif

// ************************************************************************* //
