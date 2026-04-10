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
#include "fileName.H"

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
 * \brief Export DGCompressibleFoam modal volume fields as native VTK Lagrange cells
 *
 * \details
 * dgFoamVolumeToVTU exports one native VTK high-order cell per DG cell. At
 * each write event it:
 * - looks up the conservative DG fields and thermo object from the region
 *   registry,
 * - reconstructs conservative and primitive variables directly at the native
 *   VTK Lagrange nodes,
 * - writes one VTU piece per rank and a master PVTU aggregator in parallel.
 *
 * The function object intentionally mirrors the standalone dgToVTKLagrange
 * algorithm, but uses the already-constructed solver thermo object so the
 * primitive reconstruction stays consistent with DGCompressibleFoam.
 *
 * Example declaration in \c system/controlDict:
 * \verbatim
 * functions
 * {
 *     volumeToVTU
 *     {
 *         type            dgFoamVolumeToVTU;
 *         libs            ("libDGCompressibleDerivedFunctionObjects.so");
 *         writeControl    writeTime;
 *
 *         // Optional overrides
 *         region          region0;
 *         rho             rho;
 *         rhoU            rhoU;
 *         E               E;
 *         thermo          dgThermoConservative;
 *     }
 * }
 * \endverbatim
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

    /** Guard to avoid repeating lightweight validation work */
    bool initialized_;

    /**
     * \brief Lazily validate solver-owned thermo access
     *
     * This routine only checks that the function object can resolve the
     * required mesh and thermo objects from the region registry.
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
