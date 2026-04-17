/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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
    Foam::dgMassProjection

Description
    Local DG projection driver for auxiliary spatial variables.

SourceFiles
    dgMassProjection.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgMassProjection_H
#define Foam_dgMassProjection_H

#include "List.H"
#include "Ostream.H"
#include "dgField.H"
#include "dgGeomMesh.H"
#include "error.H"
#include "pTraits.H"
#include "tensor.H"
#include "word.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class dgMassProjection Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgMassProjection
 * \brief Solve cell-local DG mass projections for auxiliary fields.
 *
 * \details
 * This class mirrors the residual-buffer part of dgTimeDiscretization, but
 * deliberately has no stage loop, saved old state, or time increment. It is
 * meant for auxiliary equations of the form
 *
 * \f[
 *     \mathbf{M}\,\hat{\mathbf{S}} = \mathbf{R},
 * \f]
 *
 * where \f$\mathbf{R}\f$ has already been assembled by a spatial DG term. The
 * class divides each residual coefficient by the cell mass-matrix diagonal,
 * writes the result into the registered field DoFs, and refreshes the matching
 * Gauss-point representation.
 *
 * Supported field/residual value types are scalar, vector, and tensor.
 */
class dgMassProjection
{
private:

    /** Helper used to trigger a clean compile-time error for unsupported types. */
    template<class Type>
    struct alwaysFalse : std::false_type
    {};

    /** Shared DG mesh used to size residual buffers and access mass matrices. */
    const dgGeomMesh& mesh_;

    /** Registered scalar auxiliary fields. */
    List<word> scalarNames_;
    List<dgField<scalar>*> scalarFields_;
    List<List<List<scalar>>> scalarResiduals_;

    /** Registered vector auxiliary fields. */
    List<word> vectorNames_;
    List<dgField<vector>*> vectorFields_;
    List<List<List<vector>>> vectorResiduals_;

    /** Registered tensor auxiliary fields. */
    List<word> tensorNames_;
    List<dgField<tensor>*> tensorFields_;
    List<List<List<tensor>>> tensorResiduals_;

    /** Return the index of \p fieldName in a name list, or -1 if absent. */
    static label findNameIndex
    (
        const List<word>& names,
        const word& fieldName
    );

    /** Validate a local DG cell index. */
    void validateCellID(const label cellID) const;

    /** Compile-time validation for supported field value types. */
    template<class Type>
    inline void checkType() const;

    /** Return the registered names for a value type. */
    template<class Type>
    inline auto& names();

    /** Return the registered names for a value type. */
    template<class Type>
    inline const auto& names() const;

    /** Return the registered field pointers for a value type. */
    template<class Type>
    inline auto& fields();

    /** Return the registered field pointers for a value type. */
    template<class Type>
    inline const auto& fields() const;

    /** Return the residual buffers for a value type. */
    template<class Type>
    inline auto& residuals();

    /** Return the residual buffers for a value type. */
    template<class Type>
    inline const auto& residuals() const;

    /** Return the residual buffer for one registered field. */
    template<class Type>
    inline List<List<Type>>& residualBuffer(const word& fieldName);

    /** Return the residual buffer for one registered field. */
    template<class Type>
    inline const List<List<Type>>& residualBuffer(const word& fieldName) const;

    /** Allocate residual storage on the mesh DoF layout. */
    template<class Type>
    void allocateResidual(List<List<Type>>& R) const;

    /** Clear one cell of one residual buffer. */
    template<class Type>
    void clearResidualCell(List<List<Type>>& R, const label cellID) const;

    /** Clear all cells of one residual buffer. */
    template<class Type>
    void clearResidual(List<List<Type>>& R) const;

    /** Solve one registered field on one cell. */
    template<class Type>
    void solveFieldCell
    (
        dgField<Type>& field,
        const List<List<Type>>& R,
        const label cellID
    ) const;

    /** Solve one registered field on all local cells. */
    template<class Type>
    void solveField
    (
        dgField<Type>& field,
        const List<List<Type>>& R
    ) const;

    /** Solve all registered fields of one value type on one cell. */
    template<class Type>
    void solveType(const label cellID) const;

    /** Solve all registered fields of one value type on all local cells. */
    template<class Type>
    void solveType() const;

    /** Write one-line residual summaries for all fields of one value type. */
    template<class Type>
    void writeResidualType(Ostream& os, const word& typeLabel) const;

public:

    // Constructors

        /**
         * \brief Construct a spatial projection driver on a DG mesh.
         *
         * \param mesh DG mesh used by the auxiliary fields
         */
        explicit dgMassProjection(const dgGeomMesh& mesh);

        /** Disallow default construction. */
        dgMassProjection() = delete;

        /** Default destructor. */
        ~dgMassProjection() = default;

    // Access

        /** \brief Return the DG mesh associated with this driver. */
        inline const dgGeomMesh& mesh() const
        {
            return mesh_;
        }

        /** \brief Return the names of registered scalar fields. */
        inline const List<word>& scalarFieldNames() const
        {
            return scalarNames_;
        }

        /** \brief Return the names of registered vector fields. */
        inline const List<word>& vectorFieldNames() const
        {
            return vectorNames_;
        }

        /** \brief Return the names of registered tensor fields. */
        inline const List<word>& tensorFieldNames() const
        {
            return tensorNames_;
        }

    // Registration

        /**
         * \brief Register an auxiliary DG field.
         *
         * The field must own DoF storage because solve() writes projected
         * modal coefficients back into the field.
         */
        template<class Type>
        void registerField(dgField<Type>& field);

    // Field/residual access

        /** \brief Return a registered scalar field by name. */
        dgField<scalar>& scalarField(const word& fieldName);

        /** \brief Return a registered scalar field by name. */
        const dgField<scalar>& scalarField(const word& fieldName) const;

        /** \brief Return a registered vector field by name. */
        dgField<vector>& vectorField(const word& fieldName);

        /** \brief Return a registered vector field by name. */
        const dgField<vector>& vectorField(const word& fieldName) const;

        /** \brief Return a registered tensor field by name. */
        dgField<tensor>& tensorField(const word& fieldName);

        /** \brief Return a registered tensor field by name. */
        const dgField<tensor>& tensorField(const word& fieldName) const;

        /** \brief Return the scalar residual buffer for a registered field. */
        List<List<scalar>>& scalarResidual(const word& fieldName);

        /** \brief Return the scalar residual buffer for a registered field. */
        const List<List<scalar>>& scalarResidual(const word& fieldName) const;

        /** \brief Return the vector residual buffer for a registered field. */
        List<List<vector>>& vectorResidual(const word& fieldName);

        /** \brief Return the vector residual buffer for a registered field. */
        const List<List<vector>>& vectorResidual(const word& fieldName) const;

        /** \brief Return the tensor residual buffer for a registered field. */
        List<List<tensor>>& tensorResidual(const word& fieldName);

        /** \brief Return the tensor residual buffer for a registered field. */
        const List<List<tensor>>& tensorResidual(const word& fieldName) const;

        /**
         * \brief Accumulate one cell residual contribution into a registered field.
         *
         * \param fieldName Registered auxiliary field name
         * \param cellID Local DG cell index
         * \param residual Residual vector indexed by local DoF
         */
        template<class Type>
        inline void addResidual
        (
            const word& fieldName,
            const label cellID,
            const List<Type>& residual
        );

    // Driver API

        /** \brief Clear all residuals of all registered fields. */
        void clearResiduals();

        /** \brief Clear all residuals on one cell for all registered fields. */
        void clearResiduals(const label cellID);

        /** \brief Clear the complete residual buffer for one registered field. */
        void clearResiduals(const word& fieldName);

        /** \brief Clear one cell residual for one registered field. */
        void clearResiduals(const word& fieldName, const label cellID);

        /** \brief Solve all registered auxiliary fields on all local cells. */
        void solve();

        /** \brief Solve all registered auxiliary fields on one local cell. */
        void solve(const label cellID);

        /** \brief Solve one registered auxiliary field on all local cells. */
        void solve(const word& fieldName);

        /** \brief Solve one registered auxiliary field on one local cell. */
        void solve(const word& fieldName, const label cellID);

    // Diagnostics

        /** \brief Write a short summary of registered auxiliary fields. */
        void writeInfo(Ostream& os) const;

        /** \brief Write maximum absolute residual magnitudes per field. */
        void writeResiduals(Ostream& os) const;
};

} // End namespace Foam

#include "dgMassProjectionI.H"

#endif

// ************************************************************************* //
