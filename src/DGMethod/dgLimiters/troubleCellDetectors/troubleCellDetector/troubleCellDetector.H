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
    Foam::troubleCellDetector

Description
    Runtime-selectable wrapper for troubled-cell detection.

SourceFiles
    troubleCellDetector.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_troubleCellDetector_H
#define Foam_troubleCellDetector_H

#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "dgField.H"
#include "dgGeomMesh.H"
#include "volFieldsFwd.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class troubleCellDetector Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::troubleCellDetector
 * \brief Runtime-selectable base class for troubled-cell sensors.
 *
 * \details
 * The wrapper stores the shared detector dictionary, the DG mesh,
 * the resolved list of fields to be checked, and optional reporting state.
 * Concrete detectors remain responsible only for any extra mandatory fields
 * that are not part of \c checkFields.
 */
class troubleCellDetector
{
public:

        /** \brief Runtime wrapper for one field listed in \c checkFields. */
        struct checkedField
        {
            /** Entry name as written in the detector dictionary. */
            word entryName;

            /** Resolved object name after applying lookupFieldName(). */
            word fieldName;

            /** True when the resolved field is a vector DG field. */
            bool isVector;

            /** Pointer to scalar DG field when \c isVector is false. */
            const dgField<scalar>* scalarFieldPtr;

            /** Pointer to vector DG field when \c isVector is true. */
            const dgField<vector>* vectorFieldPtr;
        };

protected:

        /** Local detector dictionary copied from the caller. */
        dictionary dict_;

        /** Shared DG geometric mesh used for connectivity and basis access. */
        const dgGeomMesh& mesh_;

        /** Flag enabling troubled-cell reporting to a volume indicator field. */
        bool report_;

        /** Resolved list of fields that this detector should check. */
        List<checkedField> checkFields_;

        /** Optional cell-wise 0/1 field marking troubled cells. */
        mutable autoPtr<VolumeField<label>> limitingIndicatorPtr_;

        /**
         * \brief Read an optional field-name override from the detector dictionary.
         *
         * \param dict Dictionary that may contain a custom field-name entry
         * \param entryName Name of the dictionary entry to search for
         * \param defaultName Fallback field name used when the entry is absent
         * \return The chosen field name
         */
        static word lookupFieldName
        (
            const dictionary& dict,
            const word& entryName,
            const word& defaultName
        );

        /** \brief Build the resolved check-field list from the detector dictionary. */
        void initializeCheckFields();

        /** \brief Set one cell entry in the limiting indicator if reporting is enabled. */
        void setLimitingIndicator(const label cellID, const bool flagged) const;

   public:

    /** \brief Runtime type information for the detector base class. */
    TypeName("troubleCellDetector");

    /** \brief Runtime-selection table keyed by detector type name. */
    declareRunTimeSelectionTable
    (
        autoPtr,
        troubleCellDetector,
        dictionary,
        (const dictionary& dict, const dgGeomMesh& mesh),
        (dict, mesh)
    );

    /**
     * \brief Construct the base detector and initialize shared reporting state.
     *
     * The following dictionary entries are recognized:
     * \verbatim
     * checkFields (<field1> <field2> ...);
     * report false;
     * \endverbatim
     *
     * The \c checkFields entry is required and must not be empty.
     *
     * \param dict Detector-specific dictionary
     * \param mesh DG mesh shared by the detector hierarchy
     */
    troubleCellDetector
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** \brief Disallow default construction. */
    troubleCellDetector() = delete;

    /** \brief Default virtual destructor. */
    virtual ~troubleCellDetector() = default;

    /**
     * \brief Runtime-select and construct a concrete troubled-cell detector.
     *
     * \param detectorType Name of the concrete detector type
     * \param dict Detector dictionary passed to the selected constructor
     * \param mesh DG mesh providing geometry and registry access
     * \return Auto pointer to the concrete detector
     */
    static autoPtr<troubleCellDetector> New
    (
        const word& detectorType,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** \brief Return the detector dictionary. */
    inline const dictionary& dict() const
    {
        return dict_;
    }

    /** \brief Return the DG mesh shared by this detector. */
    inline const dgGeomMesh& mesh() const
    {
        return mesh_;
    }

    /** \brief Return true when reporting to the limiting indicator is enabled. */
    inline bool report() const
    {
        return report_;
    }

    /** \brief Return the resolved list of fields to be checked. */
    inline const List<checkedField>& checkFields() const
    {
        return checkFields_;
    }

    /** \brief Reset the limiting-indicator field to zero on all cells. */
    void resetLimitingIndicator() const;

    /**
     * \brief Detect whether a cell should be marked as troubled.
     *
     * Derived classes are expected to gather any additional field, Gauss-point,
     * or modal information they need from their own cached data.
     *
     * \param cellID Target cell index
     * \return True if the cell is detected as troubled
     */
    virtual bool detect(const label cellID) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
