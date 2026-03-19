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

    The detector stores references to the conservative DG fields required by
    the current 2D detector logic:
      - Persson-Peraire uses density only
      - modified minmod uses all conservative variables

    The actual detection API receives only a cell id and derived classes
    obtain the needed GaussField data from these stored field references.

SourceFiles
    troubleCellDetector.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_troubleCellDetector_H
#define Foam_troubleCellDetector_H

#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "dgGeomMesh.H"
#include "dgField.H"

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
 * This wrapper mirrors the design used elsewhere in the codebase, such as
 * \c dgTimeScheme and \c dgThermoConservative:
 * - the concrete detector type is selected at run time from a dictionary,
 * - the base class stores the shared configuration dictionary and mesh,
 * - the base class also performs all field lookup once in the constructor.
 *
 * The field-lookup step is important because the legacy 2D detectors are not
 * purely local-in-cell. Even though the public API is simply
 * \c detect(cellID), the implementation may need to inspect:
 * - the density field \c rho for Persson-Peraire,
 * - the full conservative state \c rho, \c rhoU, and \c E for detectors
 *   derived from the modified-minmod/p-adaptive logic.
 *
 * Storing references centrally avoids repeated registry lookups inside the
 * hot detection loop and gives all derived classes a consistent way to access
 * the current Gauss-point representation through \c rhoGauss(),
 * \c rhoUGauss(), and \c EGauss().
 */
class troubleCellDetector
{
protected:

        /** Local detector dictionary copied from the caller. */
        dictionary dict_;

        /** Shared DG geometric mesh used for connectivity and basis access. */
        const dgGeomMesh& mesh_;

        /** Density field looked up from the fvMesh object registry. */
        const dgField<scalar>& rho_;

        /** Momentum field looked up from the fvMesh object registry. */
        const dgField<vector>& rhoU_;

        /** Total-energy field looked up from the fvMesh object registry. */
        const dgField<scalar>& E_;

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
     * \brief Construct the base detector and cache all required field references.
     *
     * The constructor follows the same lookup pattern used by
     * \c dgThermoConservative: fields are fetched from the fvMesh registry once
     * and then reused by every subsequent \c detect() call.
     *
     * The following optional dictionary entries are recognized:
     * \verbatim
     * rho   rho;
     * rhoU  rhoU;
     * E     E;
     * \endverbatim
     *
     * If omitted, the standard conservative field names are used.
     *
     * \param dict Detector-specific dictionary
     * \param mesh DG mesh that owns the registered DG fields
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

    /** \brief Return the cached density field reference. */
    inline const dgField<scalar>& rho() const
    {
        return rho_;
    }

    /** \brief Return the cached momentum field reference. */
    inline const dgField<vector>& rhoU() const
    {
        return rhoU_;
    }

    /** \brief Return the cached total-energy field reference. */
    inline const dgField<scalar>& E() const
    {
        return E_;
    }

    /**
     * \brief Return density Gauss data for one cell.
     *
     * \param cellID Target cell index
     * \return Cell-local GaussField view of \c rho
     */
    inline const GaussField<scalar>& rhoGauss(const label cellID) const
    {
        return rho_.gaussFields()[cellID];
    }

    /**
     * \brief Return momentum Gauss data for one cell.
     *
     * \param cellID Target cell index
     * \return Cell-local GaussField view of \c rhoU
     */
    inline const GaussField<vector>& rhoUGauss(const label cellID) const
    {
        return rhoU_.gaussFields()[cellID];
    }

    /**
     * \brief Return total-energy Gauss data for one cell.
     *
     * \param cellID Target cell index
     * \return Cell-local GaussField view of \c E
     */
    inline const GaussField<scalar>& EGauss(const label cellID) const
    {
        return E_.gaussFields()[cellID];
    }

    /**
     * \brief Detect whether a cell should be marked as troubled.
     *
     * Derived classes are expected to gather any additional Gauss-point or
     * modal information they need from the cached field references above.
     *
     * \param cellID Target cell index
     * \return True if the cell is detected as troubled
     */
    virtual bool detect(const label cellID) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
