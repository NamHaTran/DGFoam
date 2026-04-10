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

#include <vector>

namespace Foam
{
namespace
{

/**
 * \brief Point-wise conservative and derived quantities written to VTK.
 *
 * \details
 * The export utility reconstructs DG modal fields at every VTK Lagrange node.
 * For each reconstructed point it stores:
 * - conservative variables \f$\rho\f$, \f$\rho \mathbf{U}\f$, \f$E\f$,
 * - primitive variables \f$\mathbf{U}\f$, \f$p\f$, \f$T\f$,
 * - and auxiliary thermodynamic data \f$e\f$, \f$a\f$.
 *
 * The arrays are accumulated in the same order as `grid.points`, so they can
 * be passed directly to `dgVtkLagrange::writeVtuPiece()`.
 */
struct ExportPointData
{
    std::vector<scalar> rho;
    std::vector<vector> rhoU;
    std::vector<scalar> E;
    std::vector<vector> U;
    std::vector<scalar> e;
    std::vector<scalar> p;
    std::vector<scalar> T;
    std::vector<scalar> a;
};

/**
 * \brief Reconstruct a modal DG field at one interpolation point in one cell.
 *
 * \tparam Type Scalar or vector field value type.
 * \param field DG field containing modal coefficients on every cell.
 * \param cellI Cell index whose local modal expansion is evaluated.
 * \param basis Values of the local basis functions at the target point.
 * \return Reconstructed physical value \f$\sum_i \phi_i(\xi)\hat{u}_i\f$.
 *
 * \details
 * The DG solution is stored as modal coefficients per cell. This helper
 * evaluates the polynomial inside `cellI` by forming the dot product between:
 * - the modal coefficients in `field.dof()[cellI]`, and
 * - the basis values already evaluated at the local export point.
 *
 * A strict size check is performed before reconstruction because a mismatch
 * would indicate inconsistent polynomial order or an invalid nodal basis.
 */
template<class Type>
Type evaluateModalField
(
    const dgField<Type>& field,
    const label cellI,
    const List<scalar>& basis
)
{
    const cellDof<Type>& cellModes = field.dof()[cellI];

    if (cellModes.nDof() != basis.size())
    {
        FatalErrorInFunction
            << "Basis size mismatch while evaluating field " << field.name()
            << " in cell " << cellI
            << ". DoFs = " << cellModes.nDof()
            << ", basis size = " << basis.size()
            << exit(FatalError);
    }

    Type value = pTraits<Type>::zero;

    for (label modeI = 0; modeI < cellModes.nDof(); ++modeI)
    {
        value += basis[modeI]*cellModes[modeI];
    }

    return value;
}


/**
 * \brief Return the output directory for the current export time.
 */
fileName outputDir(const Time& runTime)
{
    return runTime.globalPath()/"VTK_Lagrange"/runTime.timeName();
}


/**
 * \brief Return the VTU piece file written by the current MPI rank.
 */
fileName pieceFileName(const Time& runTime)
{
    const fileName dir = outputDir(runTime);

    if (Pstream::parRun())
    {
        return dir/("dgToVTKLagrange_" + Foam::name(Pstream::myProcNo()) + ".vtu");
    }

    return dir/"dgToVTKLagrange.vtu";
}


/**
 * \brief Return the parallel PVTU manifest file for the current time.
 */
fileName pvtuFileName(const Time& runTime)
{
    return outputDir(runTime)/"dgToVTKLagrange.pvtu";
}


/**
 * \brief Create the output directory for the current time if needed.
 */
void ensureOutputDir(const Time& runTime)
{
    mkDir(outputDir(runTime));
}


/**
 * \brief Export one DG time-step to native VTK Lagrange unstructured data.
 *
 * \param runTime Current OpenFOAM time object used for naming the output.
 * \param dgMesh DG geometric mesh that defines polynomial order and cell maps.
 * \param rhoField Modal density field.
 * \param rhoUField Modal momentum field.
 * \param EField Modal total-energy field.
 * \param thermo Thermodynamic model used to recover primitive variables.
 *
 * \details
 * This routine converts a modal DG representation into a nodal VTK Lagrange
 * representation. The algorithm is:
 * 1. Determine the VTK export order corresponding to the DG polynomial order.
 * 2. Loop over all DG cells and skip null entries.
 * 3. Reject unsupported element types, while counting pyramids that are not
 *    yet exported in native VTK Lagrange form.
 * 4. For the current cell, obtain:
 *    - the PyFR standard-element points for the chosen high-order output,
 *    - the VTK node permutation (`nodemap`) expected by native Lagrange cells,
 *    - and the physical cell vertices.
 * 5. Traverse the export nodes in VTK order. For each node:
 *    - map the PyFR reference coordinate \f$\xi\f$ to the internal DG
 *      reference coordinate \f$\eta\f$,
 *    - evaluate the DG basis at \f$\eta\f$,
 *    - map \f$\eta\f$ to the physical point \f$\mathbf{x}\f$,
 *    - reconstruct \f$\rho\f$, \f$\rho\mathbf{U}\f$, and \f$E\f$ from modal
 *      coefficients,
 *    - derive \f$\mathbf{U}\f$, internal energy \f$e\f$, pressure, temperature,
 *      and speed of sound, with density clipped by `SMALL` to avoid division
 *      by zero,
 *    - append both coordinates and point-data values to contiguous buffers.
 * 6. Append cell connectivity, offsets, and the native VTK cell type so the
 *    buffered points are interpreted as one high-order Lagrange cell.
 * 7. After all cells are processed, write the VTU piece containing geometry
 *    plus scalar/vector point-data arrays.
 *
 * The resulting VTU stores solution samples at interpolation nodes, which
 * allows ParaView to render the element curvature and the high-order field
 * variation directly without first projecting to a low-order mesh.
 */
void exportTimeStep
(
    const Time& runTime,
    const dgGeomMesh& dgMesh,
    const dgField<scalar>& rhoField,
    const dgField<vector>& rhoUField,
    const dgField<scalar>& EField,
    const dgThermoConservative& thermo
)
{
    dgVtkLagrange::GridBuffers grid;
    ExportPointData pointData;
    label skippedPyramids = 0;

    const label pOrder = dgMesh.pOrder();
    const label hoOrder = dgVtkLagrange::exportOrder(pOrder);

    if (hoOrder > 8)
    {
        FatalErrorInFunction
            << "Current dgToVTKLagrange implementation supports "
            << "high-order export up to order 8. "
            << "Detected DG pOrder = " << pOrder
            << ", which maps to VTK order " << hoOrder
            << exit(FatalError);
    }

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        const dgGeomCell* cellPtr = dgMesh.cells()[cellI];

        if (!cellPtr)
        {
            continue;
        }

        const dgCellType type = cellPtr->type();

        if (!dgVtkLagrange::supportsNativeLagrange(type))
        {
            if (type == dgCellType::PYRAMID)
            {
                ++skippedPyramids;
                continue;
            }

            FatalErrorInFunction
                << "Unsupported DG cell type "
                << dgVtkLagrange::cellTypeName(type)
                << " in cell " << cellI << exit(FatalError);
        }

        const pointField& cellVertices = cellPtr->points();
        const std::vector<vector> pyfrPts =
            dgVtkLagrange::pyfrStdElementPoints(type, hoOrder);
        const std::vector<label> nodemap =
            dgVtkLagrange::vtkNodemap(type, pyfrPts.size());
        const label pointOffset = grid.points.size();

        for (const label pyfrPointI : nodemap)
        {
            const vector xi = pyfrPts[pyfrPointI];
            const vector eta = dgVtkLagrange::xiToEta(type, xi);

            List<scalar> basis;
            dgVtkLagrange::computeBasisAt(type, eta, pOrder, basis);

            const point x = mapEtaToX(eta, cellVertices, type);
            const scalar rhoVal = evaluateModalField(rhoField, cellI, basis);
            const vector rhoUVal = evaluateModalField(rhoUField, cellI, basis);
            const scalar EVal = evaluateModalField(EField, cellI, basis);

            const scalar rhoSafe = max(rhoVal, SMALL);
            const vector UVal = rhoUVal/rhoSafe;
            const scalar eVal = EVal/rhoSafe - 0.5*magSqr(UVal);

            grid.points.push_back(x);
            pointData.rho.push_back(rhoVal);
            pointData.rhoU.push_back(rhoUVal);
            pointData.E.push_back(EVal);
            pointData.U.push_back(UVal);
            pointData.e.push_back(eVal);
            pointData.T.push_back(thermo.calcTemperatureFromRhoHe(cellI, rhoSafe, eVal));
            pointData.p.push_back(thermo.calcPressureFromRhoHe(cellI, rhoSafe, eVal));
            pointData.a.push_back(thermo.calcSpeedOfSoundFromRhoHe(cellI, rhoSafe, eVal));
        }

        for (label localPointI = 0; localPointI < label(nodemap.size()); ++localPointI)
        {
            grid.connectivity.push_back(pointOffset + localPointI);
        }

        grid.offsets.push_back(grid.connectivity.size());
        grid.cellTypes.push_back(dgVtkLagrange::vtkLagrangeType(type));
        grid.partitions.push_back(Pstream::myProcNo());
    }

    const std::vector<dgVtkLagrange::ScalarPointData> scalarArrays
    {
        {"rho", &pointData.rho},
        {"E", &pointData.E},
        {"e", &pointData.e},
        {"p", &pointData.p},
        {"T", &pointData.T},
        {"a", &pointData.a}
    };

    const std::vector<dgVtkLagrange::VectorPointData> vectorArrays
    {
        {"rhoU", &pointData.rhoU},
        {"U", &pointData.U}
    };

    ensureOutputDir(runTime);
    dgVtkLagrange::writeVtuPiece
    (
        pieceFileName(runTime),
        runTime.value(),
        grid,
        scalarArrays,
        vectorArrays
    );

    if (skippedPyramids > 0)
    {
        WarningInFunction
            << "Skipped " << skippedPyramids
            << " pyramid cells because native VTK_LAGRANGE pyramid export is "
            << "not implemented yet." << nl;
    }
}

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

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

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

            exportTimeStep(runTime, dgMesh, rho, rhoU, E, thermo());
        }

        if (Pstream::parRun())
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

        Info<< "  Wrote "
            << (Pstream::parRun() ? pvtuFileName(runTime) : pieceFileName(runTime))
            << nl << endl;
    }

    return 0;
}

// ************************************************************************* //
