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
    Foam::dgFunctionObjects::dgFoamSurfaceToVTU

Description
    Export DGCompressibleFoam minus-trace values on selected boundary
    patches as VTU/PVTU surface datasets.

SourceFiles
    dgFoamSurfaceToVTU.C

\*---------------------------------------------------------------------------*/

#ifndef dgFoamSurfaceToVTU_H
#define dgFoamSurfaceToVTU_H

#include "functionObject.H"
#include "Time.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "label.H"
#include "labelList.H"
#include "scalar.H"
#include "vector.H"
#include "word.H"
#include "wordList.H"
#include "List.H"
#include "fileName.H"

namespace Foam
{

template<class Type>
class dgField;

class dgThermoConservative;

namespace dgFunctionObjects
{

/*---------------------------------------------------------------------------*\
                     Class dgFoamSurfaceToVTU Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgFunctionObjects::dgFoamSurfaceToVTU
 * \brief Export minus-trace DG values on selected patches as VTU/PVTU
 *
 * \details
 * dgFoamSurfaceToVTU writes a lightweight boundary-surface dataset for quick
 * inspection in ParaView without loading the full tessellated volume output.
 * For each selected boundary patch it:
 * - builds an equi-spaced visualization lattice on every boundary face,
 * - maps those physical surface points back to the owner-cell reference space,
 * - evaluates the conservative DG state from owner-side modal DoFs,
 * - reconstructs primitive variables from that owner-side state,
 * - writes one VTU surface piece per rank and a master PVTU file in
 *   parallel runs.
 *
 * The exported values correspond to the owner-side minus trace on the patch.
 * This makes the output inexpensive and deterministic, while avoiding any
 * dependence on boundary-condition plus states or flux-side post-processing.
 */
class dgFoamSurfaceToVTU
:
    public functionObject
{
    /** OpenFOAM time database that owns the function object */
    const Time& runTime_;

    /** Region name used to locate the target \c fvMesh */
    word regionName_;

    /** Registry name of the conservative density field */
    word rhoName_;

    /** Registry name of the conservative momentum field */
    word rhoUName_;

    /** Registry name of the conservative total-energy field */
    word EName_;

    /** Registry name of the solver thermo object */
    word thermoName_;

    /** Number of equi-spaced sub-intervals used along each face direction */
    label n_;

    /** DG polynomial order read from \c system/dgSchemes */
    label pOrder_;

    /** Merge coincident surface points to reduce duplicated patch-edge nodes */
    bool mergeSharedPoints_;

    /** Absolute tolerance used when merging coincident surface points */
    scalar pointMergeTol_;

    /** Patch names requested by the user in the function-object dictionary */
    wordList patchNames_;

    /** Resolved boundary-patch indices corresponding to \c patchNames_ */
    labelList patchIDs_;

    /** Guard to avoid repeating dictionary parsing and validation work */
    bool initialized_;

    /**
     * \brief Lazily initialize thermo constants and resolve selected patches
     *
     * This routine validates the supported thermo chain, reads the DG order
     * from \c system/dgSchemes, and converts user-provided patch names into
     * boundary-patch indices.
     */
    void initialize();

    /** Return the mesh region on which the DG solver fields are registered */
    const fvMesh& mesh() const;

    /** Return the conservative density field */
    const dgField<scalar>& rhoField() const;

    /** Return the conservative momentum field */
    const dgField<vector>& rhoUField() const;

    /** Return the conservative total-energy field */
    const dgField<scalar>& EField() const;

    /** Return the thermo object associated with the active solver region */
    const dgThermoConservative& thermo() const;

    /** Return the output directory for the current write time */
    fileName outputDir() const;

    /** Return the rank-local VTU piece path */
    fileName pieceFileName() const;

    /** Return the master PVTU aggregator path */
    fileName pvtuFileName() const;

    /**
     * \brief Evaluate a modal DG field at one sampled surface point
     *
     * \tparam Type Scalar or vector field type
     * \param field DG field to be reconstructed
     * \param cellI Owner-cell index
     * \param basis Basis-function values at the target point
     * \return Reconstructed owner-side minus-trace value
     */
    template<class Type>
    static Type evaluateField
    (
        const dgField<Type>& field,
        const label cellI,
        const List<scalar>& basis,
        const scalar highOrderScale = 1
    );

    /** Write the rank-local VTU surface piece file */
    void writePiece(const fileName& filePath) const;

    /** Write the master PVTU file that references all VTU rank pieces */
    void writePvtu(const fileName& filePath) const;

public:

    TypeName("dgFoamSurfaceToVTU");

    /**
     * \brief Construct from an OpenFOAM function-object dictionary
     *
     * The constructor stores lightweight state and immediately forwards to
     * \c read(dict) to follow the standard function-object lifecycle.
     */
    dgFoamSurfaceToVTU
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~dgFoamSurfaceToVTU() = default;

    /** Read or reload user configuration from the function-object dictionary */
    virtual bool read(const dictionary& dict) override;

    /** No-op execute hook; all export work is done in \c write() */
    virtual bool execute() override;

    /** Sample patch minus traces and write VTU/PVTU surface output */
    virtual bool write() override;
};

} // End namespace dgFunctionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
