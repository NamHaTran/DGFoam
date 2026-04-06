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
    Foam::dgLimiter

Description
    Runtime-selectable wrapper for DG limiters.

SourceFiles
    dgLimiter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgLimiter_H
#define Foam_dgLimiter_H

#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "dgField.H"
#include "dgGeomMesh.H"
#include "troubleCellDetector.H"
#include "DynamicList.H"
#include "volFieldsFwd.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgLimiter Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgLimiter
 * \brief Runtime-selectable base class for DG limiters
 *
 * \details
 * This class provides the common infrastructure used by all DG limiters in
 * the limiter pipeline.
 *
 * The base wrapper is responsible for:
 * - storing the limiter dictionary,
 * - resolving the DG fields listed in \c fields,
 * - constructing an optional troubled-cell detector,
 * - tracking which cells were limited during one pass,
 * - rebuilding cell and neighbour Gauss values after modal changes.
 *
 * The default limiter workflow is:
 * \f[
 * \text{detect troubled cells}
 * \;\rightarrow\;
 * \text{limit selected fields}
 * \;\rightarrow\;
 * \text{reconstruct Gauss values}
 * \f]
 *
 * Derived classes implement the actual limiting action through
 * \c limitField() for scalar and vector DG fields, and may override
 * \c correct() when they need to operate on every cell directly instead of
 * using a detector-driven pass.
 */
class dgLimiter
{
public:

        /** \brief Runtime wrapper for one field listed in \c fields. */
        struct limitedField
        {
            /** Resolved DG field object name. */
            word fieldName;

            /** True when the resolved field is a vector DG field. */
            bool isVector;

            /** Pointer to scalar DG field when \c isVector is false. */
            dgField<scalar>* scalarFieldPtr;

            /** Pointer to vector DG field when \c isVector is true. */
            dgField<vector>* vectorFieldPtr;
        };

protected:

        /** Resolved list of DG fields that this limiter modifies. */
        List<limitedField> limitedFields_;

        /** Runtime-selected troubled-cell detector, optional for custom limiters. */
        autoPtr<troubleCellDetector> detector_;

        /** Local limiter dictionary copied from the caller. */
        dictionary dict_;

        /** Shared DG geometric mesh used for topology and field registry access. */
        const dgGeomMesh& mesh_;

        /** Number of cells limited in the most recent correct() call. */
        label nLimitedCells_;

        /** Cells limited in the current correct() pass, in visitation order. */
        DynamicList<label> limitedCellIDs_;

        /** Flag enabling limited-cell reporting to a volume indicator field. */
        bool report_;

        /** Optional cell-wise 0/1 scalar field marking cells limited by this limiter. */
        mutable autoPtr<volScalarField> limitedCellFieldPtr_;

        /** \brief Build the resolved field list from the \c fields entry. */
        void initializeLimitedFields();

        /** \brief Construct the troubled-cell detector from the limiter dictionary when configured. */
        void initializeDetector();

        /** \brief Dispatch limiting for one resolved field in one cell. */
        void limitCell(const limitedField& fieldInfo, const label cellID);

        /** \brief Hook executed once before looping over cells in correct(). */
        virtual void preCorrect();

        /** \brief Hook executed once after looping over cells in correct(). */
        virtual void postCorrect();

        /** \brief Cache one limited cell for deferred Gauss reconstruction. */
        void cacheLimitedCell(const label cellID);

        /** \brief Set one cell entry in the limited-cell indicator field. */
        void setLimitedCellIndicator(const label cellID, const bool limited) const;

        /** \brief Reset the limited-cell indicator field to zero on all cells. */
        void resetLimitedCellIndicator() const;

        /** \brief Rebuild Gauss values on limited cells and their neighbours. */
        void rebuildLimitedCellGaussFields();

        /** \brief Limit one scalar DG field on one selected cell. */
        virtual void limitField
        (
            dgField<scalar>& field,
            const label cellID
        ) = 0;

        /** \brief Limit one vector DG field on one selected cell. */
        virtual void limitField
        (
            dgField<vector>& field,
            const label cellID
        ) = 0;

public:

    /** \brief Runtime type information for the limiter base class. */
    TypeName("dgLimiter");

    /** \brief Runtime-selection table keyed by limiter type name. */
    declareRunTimeSelectionTable
    (
        autoPtr,
        dgLimiter,
        dictionary,
        (const dictionary& dict, const dgGeomMesh& mesh),
        (dict, mesh)
    );

    /**
     * \brief Construct the base limiter and initialize shared state.
     *
     * The following dictionary entries are required:
     * \verbatim
     * type     <limiterType>;
     * fields   (<field1> <field2> ...);
     * report   false;
     * \endverbatim
     *
     * Most limiters also require
     * \verbatim
     * detector <detectorType>;
     * \endverbatim
     * plus a detector-coefficients sub-dictionary named
     * \c <detectorType>Coeffs, for example \c PerssonPeraireCoeffs.
     * Specialized derived limiters may override \c correct() and operate
     * without a detector.
     *
     * \param dict Limiter dictionary
     * \param mesh DG mesh shared by the limiter hierarchy
     */
    dgLimiter
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** \brief Disallow default construction. */
    dgLimiter() = delete;

    /** \brief Default virtual destructor. */
    virtual ~dgLimiter() = default;

    /**
     * \brief Runtime-select and construct a concrete limiter from \c dict.type.
     *
     * \param dict Limiter dictionary
     * \param mesh DG mesh providing geometry and registry access
     * \return Auto pointer to the concrete limiter
     */
    static autoPtr<dgLimiter> New
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** \brief Re-read the limiter dictionary and rebuild shared state. */
    virtual void read(const dictionary& dict);

    /** \brief Apply limiting using the configured trouble-cell detector. */
    virtual void correct();

    /** \brief Return the limiter dictionary. */
    inline const dictionary& dict() const
    {
        return dict_;
    }

    /** \brief Return the DG mesh shared by this limiter. */
    inline const dgGeomMesh& mesh() const
    {
        return mesh_;
    }

    /** \brief Return the resolved list of DG fields to be limited. */
    inline const List<limitedField>& limitedFields() const
    {
        return limitedFields_;
    }

    /** \brief Return true when this limiter has a configured detector. */
    inline bool hasDetector() const
    {
        return detector_.valid();
    }

    /** \brief Return the troubled-cell detector used by this limiter. */
    inline const troubleCellDetector& detector() const
    {
        return detector_();
    }

    /** \brief Return the number of cells limited in the last correct() call. */
    inline label nLimitedCells() const
    {
        return nLimitedCells_;
    }

    /** \brief Return true when limited-cell reporting is enabled. */
    inline bool report() const
    {
        return report_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
