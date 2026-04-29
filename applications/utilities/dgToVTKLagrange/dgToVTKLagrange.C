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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Export DGFoam modal DG results to native VTK Lagrange VTU/PVTU files.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "OSspecific.H"
#include "Pstream.H"

#include "dgCellType.H"
#include "dgField.H"
#include "dgAffineMapping.H"
#include "dgGeomMesh.H"
#include "dgThermoConservative.H"
#include "dgVtkLagrangeTools.H"
#include "basisFunctions.H"

#include <array>
#include <vector>

namespace Foam
{
namespace
{


#include "dgToVTKLagrangeTypes.H"
#include "dgToVTKLagrangeFieldOps.H"
#include "dgToVTKLagrangeOutput.H"
#include "dgToVTKLagrangeExportVolume.H"
#include "dgToVTKLagrangeExportPatches.H"

} // End anonymous namespace
} // End namespace Foam


using namespace Foam;

/**
 * \brief Program entry point for DG-to-VTK Lagrange export.
 *
 * \details
 * The executable selects requested time directories, rebuilds the DG fields for
 * each time, and delegates the per-time export to `exportTimeStep()`.
 *
 * In serial mode it writes a single `.vtu` file per time. In parallel mode
 * every rank writes one VTU piece and the master rank additionally writes a
 * `.pvtu` manifest that references all rank-local pieces and declares the
 * exported point-data schema.
 */
int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, true);

    argList::addNote
    (
        "Export DGFoam modal results to VTK_LAGRANGE VTU/PVTU files."
    );
    argList::addBoolOption
    (
        "no-internal",
        "Do not export volume cells; write only high-order boundary patches"
    );
    argList::addOption
    (
        "patchList",
        "wordList",
        "Boundary patches to export with -no-internal, e.g. \"(inlet wall)\""
    );
    argList::addOption
    (
        "lagrangeLimiterTolerance",
        "scalar",
        "Tolerance used by the internal Lagrange-node positivity limiter "
        "(default: 1e-6)"
    );
    argList::addOption
    (
        "smoothByClampingTheta",
        "scalar",
        "If theta1/theta2 fields exist, zero all high-order modes of "
        "rho/rhoU/E in cells where inputTheta <= theta1 < 1 or "
        "inputTheta <= theta2 < 1, and in their local neighbour cells"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const scalar lagrangeLimiterTolerance =
        args.getOrDefault<scalar>("lagrangeLimiterTolerance", 1.0e-6);
    const bool noInternal = args.found("no-internal");
    const bool hasPatchList = args.found("patchList");
    const wordList patchNames =
        hasPatchList ? args.getList<word>("patchList") : wordList();
    const bool smoothByClampingTheta = args.found("smoothByClampingTheta");
    const scalar clampTheta =
        smoothByClampingTheta
      ? args.get<scalar>("smoothByClampingTheta")
      : scalar(1);

    if (lagrangeLimiterTolerance <= 0)
    {
        FatalErrorInFunction
            << "Option -lagrangeLimiterTolerance must be positive, but got "
            << lagrangeLimiterTolerance << exit(FatalError);
    }

    if (smoothByClampingTheta && (clampTheta < 0 || clampTheta > 1))
    {
        FatalErrorInFunction
            << "Option -smoothByClampingTheta must be in [0, 1], but got "
            << clampTheta << exit(FatalError);
    }

    if (hasPatchList && !noInternal)
    {
        FatalErrorInFunction
            << "Option -patchList is only valid together with -no-internal."
            << exit(FatalError);
    }

    dgGeomMesh dgMesh(mesh);

    Info<< "Create DG geometric mesh with polynomial order "
        << dgMesh.pOrder() << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Exporting time " << runTime.timeName() << nl;

        {
            dgField<scalar> rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                true
            );

            dgField<vector> rhoU
            (
                IOobject
                (
                    "rhoU",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                true
            );

            dgField<scalar> E
            (
                IOobject
                (
                    "E",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                true
            );

            autoPtr<volScalarField> theta1Ptr(nullptr);
            autoPtr<volScalarField> theta2Ptr(nullptr);
            bool activeSmoothByClampingTheta = smoothByClampingTheta;

            if (smoothByClampingTheta)
            {
                IOobject theta1IO
                (
                    "theta1",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject theta2IO
                (
                    "theta2",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if
                (
                    theta1IO.typeHeaderOk<volScalarField>(false)
                 && theta2IO.typeHeaderOk<volScalarField>(false)
                )
                {
                    theta1Ptr.reset(new volScalarField(theta1IO, mesh));
                    theta2Ptr.reset(new volScalarField(theta2IO, mesh));
                }
                else
                {
                    activeSmoothByClampingTheta = false;

                    WarningInFunction
                        << "Option -smoothByClampingTheta was requested, but "
                        << "could not find both positivity-preserving limiter "
                        << "fields theta1 and theta2 at time "
                        << runTime.timeName()
                        << ". Disabling smoothByClampingTheta for this time "
                        << "and exporting normally." << nl;
                }
            }

            dgField<scalar> p
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                false
            );

            dgField<scalar> T
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                false
            );

            dgField<vector> U
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                false
            );

            dgField<vector> SRho
            (
                IOobject
                (
                    "SRho",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                true
            );

            dgField<tensor> SRhoU
            (
                IOobject
                (
                    "SRhoU",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                true
            );

            dgField<vector> SE
            (
                IOobject
                (
                    "SE",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                true
            );

            dgField<tensor> gradU
            (
                IOobject
                (
                    "gradU",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                false
            );

            dgField<vector> gradP
            (
                IOobject
                (
                    "gradP",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                false
            );

            dgField<vector> gradT
            (
                IOobject
                (
                    "gradT",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dgMesh,
                false
            );

            IOdictionary thermoDict
            (
                IOobject
                (
                    "thermophysicalProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            const dictionary& dgThermoDict = thermoDict.subDict("dgThermo");
            const word dgThermoTypeName = dgThermoDict.get<word>("type");

            autoPtr<dgThermoConservative> thermo =
                dgThermoConservative::New(dgThermoTypeName, thermoDict, dgMesh);

            if (noInternal)
            {
                exportPatchOnlyTimeStep
                (
                    runTime,
                    mesh,
                    dgMesh,
                    rho,
                    rhoU,
                    E,
                    thermo(),
                    lagrangeLimiterTolerance,
                    activeSmoothByClampingTheta,
                    clampTheta,
                    theta1Ptr.valid() ? &theta1Ptr() : nullptr,
                    theta2Ptr.valid() ? &theta2Ptr() : nullptr,
                    patchNames
                );
            }
            else
            {
                exportTimeStep
                (
                    runTime,
                    dgMesh,
                    rho,
                    rhoU,
                    E,
                    thermo(),
                    lagrangeLimiterTolerance,
                    activeSmoothByClampingTheta,
                    clampTheta,
                    theta1Ptr.valid() ? &theta1Ptr() : nullptr,
                    theta2Ptr.valid() ? &theta2Ptr() : nullptr
                );
            }
        }

        if (!noInternal && Pstream::parRun())
        {
            Pstream::barrier(UPstream::worldComm);

            if (Pstream::master())
            {
                const std::vector<dgVtkLagrange::PointDataSchema> pointData
                {
                    {"rho", 1},
                    {"rhoU", 3},
                    {"E", 1},
                    {"U", 3},
                    {"e", 1},
                    {"p", 1},
                    {"T", 1},
                    {"a", 1}
                };

                ensureOutputDir(runTime);
                dgVtkLagrange::writePvtu
                (
                    pvtuFileName(runTime),
                    runTime.value(),
                    Pstream::nProcs(),
                    "dgToVTKLagrange",
                    pointData
                );
            }

            Pstream::barrier(UPstream::worldComm);
        }

    }

    if (!noInternal && (!Pstream::parRun() || Pstream::master()))
    {
        writePvdManifest(runTime, timeDirs);
    }

    return 0;
}

// ************************************************************************* //
