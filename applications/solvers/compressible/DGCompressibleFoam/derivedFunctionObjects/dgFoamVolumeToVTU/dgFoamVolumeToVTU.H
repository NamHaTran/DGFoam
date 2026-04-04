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
    Foam::dgFunctionObjects::dgFoamVolumeToVTU

Description
    Export DGCompressibleFoam volume samples to VTU for ParaView.

SourceFiles
    dgFoamVolumeToVTU.C

\*---------------------------------------------------------------------------*/

#ifndef dgFoamVolumeToVTU_H
#define dgFoamVolumeToVTU_H

#include "functionObject.H"
#include "Time.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "word.H"
#include "List.H"
#include "fileName.H"
#include "dgEquiSpacedTessellation.H"

#include <memory>

namespace Foam
{

template<class Type>
class dgField;

class dgThermoConservative;

namespace dgFunctionObjects
{

/*---------------------------------------------------------------------------*\
                      Class dgFoamVolumeToVTU Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgFunctionObjects::dgFoamVolumeToVTU
 * \brief Sample DGCompressibleFoam volume fields and write VTU/PVTU output
 *
 * \details
 * dgFoamVolumeToVTU exports a tessellated volumetric representation of the
 * reconstructed DG solution for ParaView. At each write event it:
 * - looks up the conservative DG fields and thermo object from the region
 *   registry,
 * - builds an equi-spaced sampling lattice in the reference element of each
 *   supported cell topology,
 * - evaluates modal DoFs at those sampling points,
 * - reconstructs primitive variables pointwise from the sampled conservative
 *   state,
 * - writes one VTU piece per rank and a master PVTU aggregator in parallel.
 *
 * The generated dataset is intended for visualization rather than restart or
 * exact numerical post-processing. Optional point merging can be used to
 * smooth element-interface jumps for faster qualitative inspection in
 * ParaView.
 */
class dgFoamVolumeToVTU
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

    /** Number of equi-spaced sub-intervals along each reference direction */
    label n_;

    /** DG polynomial order read from \c system/dgSchemes */
    label pOrder_;

    /** Enable collapsed through-thickness output for extruded 2D cases */
    bool twoDMode_;

    /** Cartesian axis normal to the 2D plane when \c twoDMode_ is enabled */
    label twoDNormalAxis_;

    /** Merge coincident tessellation points for smoother visualization */
    bool mergeSharedPoints_;

    /** Absolute tolerance used when merging coincident sampled points */
    scalar pointMergeTol_;

    /** Constant-pressure heat capacity read from thermo dictionaries */
    scalar Cp_;

    /** Constant-volume heat capacity reconstructed as \c Cp - R */
    scalar Cv_;

    /** Specific gas constant taken from the active equation of state */
    scalar R_;

    /** Guard to avoid repeating dictionary parsing and validation work */
    bool initialized_;

    /** Cached sampling/tessellation helper shared across supported cell types */
    std::unique_ptr<dgEquiSpacedTessellation> tessellationPtr_;

    /**
     * \brief Lazily initialize thermo constants and sampling configuration
     *
     * This routine validates the supported thermo chain, reads the DG order
     * from \c system/dgSchemes, and prepares the reusable tessellation helper.
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

    /** Return the cached equi-spaced tessellation helper */
    const dgEquiSpacedTessellation& tessellation() const;

    /** Return the output directory for the current write time */
    fileName outputDir() const;

    /** Return the rank-local VTU piece path */
    fileName pieceFileName() const;

    /** Return the master PVTU aggregator path */
    fileName pvtuFileName() const;

    /**
     * \brief Evaluate a modal DG field at one sampling point
     *
     * \tparam Type Scalar or vector field type
     * \param field DG field to be reconstructed
     * \param cellI Cell index
     * \param basis Basis-function values at the target point
     * \return Reconstructed field value at the sampling point
     */
    template<class Type>
    static Type evaluateField
    (
        const dgField<Type>& field,
        const label cellI,
        const List<scalar>& basis
    );

    /** Write the rank-local VTU piece file */
    void writePiece(const fileName& filePath) const;

    /** Write the master PVTU file that references all VTU rank pieces */
    void writePvtu(const fileName& filePath) const;

public:

    TypeName("dgFoamVolumeToVTU");

    /**
     * \brief Construct from an OpenFOAM function-object dictionary
     *
     * The constructor stores lightweight state and immediately forwards to
     * \c read(dict) to follow the standard function-object lifecycle.
     */
    dgFoamVolumeToVTU
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~dgFoamVolumeToVTU() = default;

    /** Read or reload user configuration from the function-object dictionary */
    virtual bool read(const dictionary& dict) override;

    /** No-op execute hook; all export work is done in \c write() */
    virtual bool execute() override;

    /** Sample DG volume fields and write VTU/PVTU output */
    virtual bool write() override;
};

} // End namespace dgFunctionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
