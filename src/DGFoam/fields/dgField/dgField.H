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

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgField_H
#define Foam_dgField_H

#include "regIOobject.H"
#include "dofField.H"
#include "GaussField.H"
#include "dgGeomMesh.H"
#include "dgGeneralBoundaryManager.H"
#include "tensor.H"

namespace Foam
{

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

        /** Reference to the DG geometric mesh */
        const dgGeomMesh& mesh_;

        /** Flag indicating whether this field owns DoF data */
        bool hasDof_;

        /** Pointer to DoF field (null for DoF-less fields) */
        autoPtr<dofField<Type>> dofPtr_;

        /** GaussField values for each cell in the mesh */
        List<GaussField<Type>> gaussFields_;

public:

    //- Runtime type information
    TypeName("dgField");

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct DG field
         *
         * Creates a DG field associated with a DG geometric mesh.
         * Optionally allocates DoF storage.
         */
        dgField
        (
            const IOobject& io,
            const dgGeomMesh& mesh,
            bool hasDof = true
        );

        /** Copy construction is disabled */
        dgField(const dgField&) = delete;

        /** Copy assignment is disabled */
        dgField& operator=(const dgField&) = delete;

    //====================================================================//
    //                           Access
    //====================================================================//

        /** Return true if this field owns DoF data */
        inline bool hasDof() const
        {
            return dofPtr_ != nullptr;
        }

        /**
         * \brief Return non-const access to DoF field
         *
         * Fatal error if this field does not own DoF data.
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
         * \brief Return const access to DoF field
         *
         * Fatal error if this field does not own DoF data.
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

        /** Return non-const access to Gauss fields */
        inline List<GaussField<Type>>& gaussFields()
        {
            return gaussFields_;
        }

        /** Return const access to Gauss fields */
        inline const List<GaussField<Type>>& gaussFields() const
        {
            return gaussFields_;
        }

        /** Return number of cells */
        inline label nCells() const
        {
            return gaussFields_.size();
        }

    //====================================================================//
    //                              I/O
    //====================================================================//

        /**
         * \brief Write field data
         *
         * Writes field information via regIOobject.
         */
        bool writeData(Ostream& os) const;

    //====================================================================//
    //                           Assignment
    //====================================================================//

        /**
         * \brief Assign uniform value to all Gauss fields
         */
        dgField<Type>& operator=(const Type& value);
};

} // End namespace Foam

#include "dgFieldI.H"

#endif

// ************************************************************************* //
