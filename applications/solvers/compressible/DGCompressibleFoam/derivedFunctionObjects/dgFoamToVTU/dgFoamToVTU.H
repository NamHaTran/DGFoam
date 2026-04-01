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
    Foam::dgFunctionObjects::dgFoamToVTU

Description
    Export DGCompressibleFoam solution samples to VTU for ParaView.

    The function object samples the conservative DG fields \c rho,
    \c rhoU, and \c E at a denser set of reference-space points inside
    every cell, reconstructs the primitive quantities pointwise, and
    writes the resulting point cloud as a VTU/PVTU dataset.

    The current implementation is intentionally solver-specific:
    - it assumes the fields and thermo object created by
      \c DGCompressibleFoam are present in the region database,
    - it reads the equi-spaced output density from dictionary entry
      \c outputPoints (or legacy alias \c n),
    - it currently supports the thermo chain
      \c idealGas + \c constantCp + \c sensibleInternalEnergy.

    The exported VTU contains point data arrays:
    - \c p : pressure
    - \c T : temperature
    - \c U : velocity

    The sampled points are also connected into a linear sub-cell
    tessellation so ParaView can render a volumetric representation of the
    reconstructed DG solution instead of only a point cloud.

    Optional smoothing for visualization can be enabled with
    \c mergeSharedPoints, which merges coincident sampled points across
    neighboring cells and averages their conservative state before
    reconstructing primitive variables. This produces more continuous-looking
    contours, at the cost of no longer showing the exact DG jump at element
    interfaces.

    For extruded 2D cases, \c twoDNormalAxis can be set to \c x, \c y, or
    \c z. In that mode the exporter collapses the sampling lattice to a
    single layer through the specified out-of-plane direction and writes
    in-plane surface cells while keeping their 3D coordinates.

\*---------------------------------------------------------------------------*/

#ifndef dgFoamToVTU_H
#define dgFoamToVTU_H

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
                           Class dgFoamToVTU Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgFunctionObjects::dgFoamToVTU
 * \brief Sample DGCompressibleFoam conservative fields and write VTU output
 *
 * \details
 * dgFoamToVTU is a lightweight \c functionObject that piggybacks on the
 * OpenFOAM runtime execution/writing controls. At each write event it:
 * - looks up the DG conservative fields from the mesh registry,
 * - constructs a cached sampling lattice in reference space for each cell
 *   topology,
 * - evaluates modal DG fields at those sample points using the cell DoFs,
 * - computes \f$\mathbf{U}\f$, \f$T\f$, and \f$p\f$ pointwise from the
 *   sampled conservative state,
 * - writes a VTU piece file for the local rank and an optional PVTU file
 *   on the master rank.
 *
 * The generated dataset is a tessellated unstructured grid:
 * - HEX cells are subdivided into linear hexahedra,
 * - PRISM cells are subdivided into linear wedges,
 * - PYRAMID cells are subdivided into a mix of hexahedra, wedges, and
 *   pyramids,
 * - TET cells are subdivided into a mix of wedges and tetrahedra.
 *
 * No external meshing library is required. The tessellation is generated
 * directly from simple reference-space connectivity templates.
 */
class dgFoamToVTU
:
    public functionObject
{
    /** Reference to the OpenFOAM time database that owns this function object */
    const Time& runTime_;

    /** Region name used to find the target \c fvMesh */
    word regionName_;

    /** Registry name of the conservative density field */
    word rhoName_;

    /** Registry name of the conservative momentum field */
    word rhoUName_;

    /** Registry name of the conservative total-energy field */
    word EName_;

    /** Registry name of the solver thermo object */
    word thermoName_;

    /**
     * Internal number of equi-spaced sub-intervals per reference direction.
     * The user-facing output point count is therefore \c n_ + 1.
     */
    label n_;

    /** DG polynomial order read from \c system/dgSchemes */
    label pOrder_;

    /** Activate collapsed through-thickness sampling for extruded 2D cases */
    bool twoDMode_;

    /** Cartesian component normal to the 2D plane: 0=x, 1=y, 2=z */
    label twoDNormalAxis_;

    /** Merge coincident tessellation points to obtain smoother contours */
    bool mergeSharedPoints_;

    /** Absolute tolerance used when merging coincident sample points */
    scalar pointMergeTol_;

    /** Constant-pressure heat capacity read from thermo dictionaries */
    scalar Cp_;

    /** Constant-volume heat capacity reconstructed as \c Cp - R */
    scalar Cv_;

    /** Specific gas constant taken from the active EOS model */
    scalar R_;

    /** Guard to ensure expensive dictionary/model validation is done once */
    bool initialized_;

    /** Generic sampling and tessellation helper shared across cell types */
    std::unique_ptr<dgEquiSpacedTessellation> tessellationPtr_;

    /**
     * \brief Lazily read supporting dictionaries and validate the thermo chain
     *
     * This routine also extracts the scalar constants needed for pointwise
     * conservative-to-primitive reconstruction.
     */
    void initialize();

    /** Return the mesh region on which the solver fields are registered */
    const fvMesh& mesh() const;

    /** Return the conservative density DG field */
    const dgField<scalar>& rhoField() const;

    /** Return the conservative momentum DG field */
    const dgField<vector>& rhoUField() const;

    /** Return the conservative total-energy DG field */
    const dgField<scalar>& EField() const;

    /** Return the thermo object associated with the active solver region */
    const dgThermoConservative& thermo() const;

    /** Return the generic equi-spaced tessellation helper */
    const dgEquiSpacedTessellation& tessellation() const;

    /** Return the post-processing directory for the current time */
    fileName outputDir() const;

    /** Return the VTU file path for the current rank */
    fileName pieceFileName() const;

    /** Return the PVTU aggregator file path written on the master rank */
    fileName pvtuFileName() const;

    /**
     * \brief Evaluate a DG field at one sample point from modal DoFs
     *
     * \tparam Type Scalar or vector field type
     * \param field DG field to be evaluated
     * \param cellI Cell index
     * \param basis Precomputed modal basis values at the target point
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

    TypeName("dgFoamToVTU");

    /**
     * \brief Construct from the function-object dictionary
     *
     * The constructor stores lightweight configuration and immediately calls
     * \c read(dict) so the object follows the normal OpenFOAM function-object
     * lifecycle.
     */
    dgFoamToVTU
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~dgFoamToVTU() = default;

    /** Read or reload user configuration from dictionary entries */
    virtual bool read(const dictionary& dict) override;

    /** No-op execute hook; the exporter performs all work in \c write() */
    virtual bool execute() override;

    /** Sample DG fields and write VTU/PVTU output for the current time */
    virtual bool write() override;
};

} // End namespace dgFunctionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
