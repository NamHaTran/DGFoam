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
#include "troubleCellDetector.H"

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


/** \brief Conservative state sampled at one Lagrange node. */
struct ConservativeSample
{
    scalar rho;
    vector rhoU;
    scalar E;
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
 * \brief Return the zeroth modal coefficient for one DG field in one cell.
 *
 * \details
 * The modal basis used by DGFoam stores the cell-average value in DoF0. The
 * \c smoothContinuousFields exporter option uses this value as a constant
 * reconstruction in cells that are not marked troubled by KXRCF.
 */
template<class Type>
Type evaluateCellMeanField
(
    const dgField<Type>& field,
    const label cellI
)
{
    const cellDof<Type>& cellModes = field.dof()[cellI];

    if (cellModes.nDof() == 0)
    {
        FatalErrorInFunction
            << "Field " << field.name() << " has no DoF in cell "
            << cellI << exit(FatalError);
    }

    return cellModes[0];
}


/**
 * \brief Zhang-Shu limiter coefficient clamped to [0, 1].
 */
scalar limiterTheta
(
    const scalar meanValue,
    const scalar minValue,
    const scalar omega,
    const scalar epsilon
)
{
    const scalar scale = max(scalar(1), mag(meanValue));

    if (mag(meanValue - minValue) <= epsilon*scale)
    {
        return 1.0;
    }

    const scalar rawTheta = (meanValue - omega)/(meanValue - minValue);
    return min(max(rawTheta, scalar(0)), scalar(1));
}


/**
 * \brief Return rho*e from conservative variables with a guarded rho.
 */
scalar thermoEnergyDensity
(
    const scalar rho,
    const vector& rhoU,
    const scalar E,
    const scalar epsilon
)
{
    const scalar rhoSafe =
        (mag(rho) > epsilon ? rho : (rho >= 0 ? epsilon : -epsilon));

    return E - 0.5*magSqr(rhoU)/rhoSafe;
}


/**
 * \brief Scale the non-mean modal content in one cell.
 */
template<class Type>
void scaleHighOrderModes
(
    dgField<Type>& field,
    const label cellI,
    const scalar theta
)
{
    if (!field.hasDof() || theta >= 1.0)
    {
        return;
    }

    cellDof<Type>& cellModes = field.dof()[cellI];

    for (label modeI = 1; modeI < cellModes.nDof(); ++modeI)
    {
        cellModes[modeI] *= theta;
    }

    field.dof().updateCellDof(cellI);
}


/**
 * \brief Limit modal conservative fields using the same Lagrange nodes that are
 *        later written to VTK.
 */
void limitConservativeFieldsAtLagrangeNodes
(
    const dgGeomMesh& dgMesh,
    dgField<scalar>& rhoField,
    dgField<vector>& rhoUField,
    dgField<scalar>& EField,
    const scalar epsilon
)
{
    label nLimitedCells = 0;

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

        if (rhoField.dof()[cellI].nDof() <= 1)
        {
            continue;
        }

        const dgCellType type = cellPtr->type();

        if (!dgVtkLagrange::supportsNativeLagrange(type))
        {
            if (type == dgCellType::PYRAMID)
            {
                continue;
            }

            FatalErrorInFunction
                << "Unsupported DG cell type "
                << dgVtkLagrange::cellTypeName(type)
                << " in cell " << cellI << exit(FatalError);
        }

        const std::vector<vector> pyfrPts =
            dgVtkLagrange::pyfrStdElementPoints(type, hoOrder);
        const std::vector<label> nodemap =
            dgVtkLagrange::vtkNodemap(type, pyfrPts.size());

        DynamicList<ConservativeSample> samples(nodemap.size());
        scalar minRho = GREAT;

        for (const label pyfrPointI : nodemap)
        {
            const vector eta =
                dgVtkLagrange::xiToEta(type, pyfrPts[pyfrPointI]);

            List<scalar> basis;
            dgVtkLagrange::computeBasisAt(type, eta, pOrder, basis);

            ConservativeSample sample;
            sample.rho = evaluateModalField(rhoField, cellI, basis);
            sample.rhoU = evaluateModalField(rhoUField, cellI, basis);
            sample.E = evaluateModalField(EField, cellI, basis);

            minRho = min(minRho, sample.rho);
            samples.append(sample);
        }

        const scalar meanRho = rhoField.dof()[cellI].dof()[0];
        const vector& meanRhoU = rhoUField.dof()[cellI].dof()[0];
        const scalar meanE = EField.dof()[cellI].dof()[0];

        const scalar theta1 =
            limiterTheta(meanRho, minRho, min(epsilon, meanRho), epsilon);

        scalar minRhoThermoEnergy = GREAT;

        forAll(samples, sampleI)
        {
            const ConservativeSample& sample = samples[sampleI];
            const scalar limitedRho =
                meanRho + theta1*(sample.rho - meanRho);

            minRhoThermoEnergy =
                min
                (
                    minRhoThermoEnergy,
                    thermoEnergyDensity
                    (
                        limitedRho,
                        sample.rhoU,
                        sample.E,
                        epsilon
                    )
                );
        }

        const scalar meanRhoThermoEnergy =
            thermoEnergyDensity(meanRho, meanRhoU, meanE, epsilon);

        const scalar theta2 =
            limiterTheta
            (
                meanRhoThermoEnergy,
                minRhoThermoEnergy,
                min(min(epsilon, meanRho), meanRhoThermoEnergy),
                epsilon
            );

        if (theta1 >= 1.0 && theta2 >= 1.0)
        {
            continue;
        }

        scaleHighOrderModes(rhoField, cellI, theta1*theta2);
        scaleHighOrderModes(rhoUField, cellI, theta2);
        scaleHighOrderModes(EField, cellI, theta2);

        ++nLimitedCells;
    }

    reduce(nLimitedCells, sumOp<label>());

    Info<< "Lagrange-node positivity limiter applied on "
        << nLimitedCells << " cell(s)." << nl;
}


/**
 * \brief Detect troubled cells using the KXRCF inflow-jump sensor.
 */
List<bool> detectKXRCFTroubledCells
(
    const dgGeomMesh& dgMesh,
    const scalar threshold,
    const bool LPRMode,
    label& nTroubledCells
)
{
    dictionary detectorDict;

    wordList checkFields(3);
    checkFields[0] = "rho";
    checkFields[1] = "rhoU";
    checkFields[2] = "E";

    detectorDict.set("checkFields", checkFields);
    detectorDict.set("directionField", word("rhoU"));
    detectorDict.set("threshold", threshold);
    detectorDict.set("LPRMode", LPRMode);
    detectorDict.set("includeBoundaryFaces", true);
    detectorDict.set("report", false);

    autoPtr<troubleCellDetector> detector =
        troubleCellDetector::New("KXRCF", detectorDict, dgMesh);

    List<bool> troubledCells(dgMesh.nCells(), false);
    nTroubledCells = 0;

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        if (!dgMesh.cells()[cellI])
        {
            continue;
        }

        if (detector->detect(cellI))
        {
            troubledCells[cellI] = true;
            ++nTroubledCells;
        }
    }

    reduce(nTroubledCells, sumOp<label>());

    return troubledCells;
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
 * 2. Optionally use the KXRCF detector to mark troubled cells when
 *    \c smoothContinuousFields is enabled.
 * 3. Apply the internal Lagrange-node positivity limiter.
 * 4. Loop over all DG cells and skip null entries.
 * 5. Reject unsupported element types, while counting pyramids that are not
 *    yet exported in native VTK Lagrange form.
 * 6. For the current cell, obtain:
 *    - the PyFR standard-element points for the chosen high-order output,
 *    - the VTK node permutation (`nodemap`) expected by native Lagrange cells,
 *    - and the physical cell vertices.
 * 7. Traverse the export nodes in VTK order. For each node:
 *    - map the PyFR reference coordinate \f$\xi\f$ to the internal DG
 *      reference coordinate \f$\eta\f$,
 *    - map \f$\eta\f$ to the physical point \f$\mathbf{x}\f$,
 *    - reconstruct \f$\rho\f$, \f$\rho\mathbf{U}\f$, and \f$E\f$ from modal
 *      coefficients,
 *    - if \c smoothContinuousFields is enabled, reconstruct with full DoF only
 *      on KXRCF-marked cells and use DoF0 as a constant value elsewhere,
 *    - derive \f$\mathbf{U}\f$, internal energy \f$e\f$, pressure, temperature,
 *      and speed of sound, with density clipped by `SMALL` to avoid division
 *      by zero,
 *    - append both coordinates and point-data values to contiguous buffers.
 * 8. Append cell connectivity, offsets, and the native VTK cell type so the
 *    buffered points are interpreted as one high-order Lagrange cell.
 * 9. After all cells are processed, write the VTU piece containing geometry
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
    dgField<scalar>& rhoField,
    dgField<vector>& rhoUField,
    dgField<scalar>& EField,
    const dgThermoConservative& thermo,
    const scalar limiterEpsilon,
    const bool smoothContinuousFields,
    const scalar kxrcfThreshold,
    const bool kxrcfLPRMode
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

    List<bool> troubledCells(dgMesh.nCells(), true);
    label nTroubledCells = 0;

    if (smoothContinuousFields)
    {
        troubledCells =
            detectKXRCFTroubledCells
            (
                dgMesh,
                kxrcfThreshold,
                kxrcfLPRMode,
                nTroubledCells
            );
    }

    Info<< "KXRCF troubled cells marked for full-DoF Lagrange "
        << "reconstruction: " << nTroubledCells;

    if (!smoothContinuousFields)
    {
        Info<< " (smoothContinuousFields disabled)";
    }

    Info<< nl;

    limitConservativeFieldsAtLagrangeNodes
    (
        dgMesh,
        rhoField,
        rhoUField,
        EField,
        limiterEpsilon
    );

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
        const bool useFullDof =
            !smoothContinuousFields || troubledCells[cellI];

        const scalar rhoMean = evaluateCellMeanField(rhoField, cellI);
        const vector rhoUMean = evaluateCellMeanField(rhoUField, cellI);
        const scalar EMean = evaluateCellMeanField(EField, cellI);

        for (const label pyfrPointI : nodemap)
        {
            const vector xi = pyfrPts[pyfrPointI];
            const vector eta = dgVtkLagrange::xiToEta(type, xi);

            const point x = mapEtaToX(eta, cellVertices, type);

            scalar rhoVal = rhoMean;
            vector rhoUVal = rhoUMean;
            scalar EVal = EMean;

            if (useFullDof)
            {
                List<scalar> basis;
                dgVtkLagrange::computeBasisAt(type, eta, pOrder, basis);

                rhoVal = evaluateModalField(rhoField, cellI, basis);
                rhoUVal = evaluateModalField(rhoUField, cellI, basis);
                EVal = evaluateModalField(EField, cellI, basis);
            }

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
    argList::addOption
    (
        "lagrangeLimiterTolerance",
        "scalar",
        "Tolerance used by the internal Lagrange-node positivity limiter "
        "(default: 1e-6)"
    );
    argList::addOption
    (
        "kxrcfThreshold",
        "scalar",
        "KXRCF troubled-cell detector threshold used by "
        "-smoothContinuousFields (default: 1)"
    );
    argList::addBoolOption
    (
        "kxrcfLPRMode",
        "Use the LPR normalization in the KXRCF detector when "
        "-smoothContinuousFields is enabled"
    );
    argList::addBoolOption
    (
        "smoothContinuousFields",
        "Use KXRCF to keep full DoF only near strong discontinuities and "
        "write DoF0 values in smooth cells"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const scalar lagrangeLimiterTolerance =
        args.getOrDefault<scalar>("lagrangeLimiterTolerance", 1.0e-6);
    const scalar kxrcfThreshold =
        args.getOrDefault<scalar>("kxrcfThreshold", 1.0);
    const bool kxrcfLPRMode = args.found("kxrcfLPRMode");
    const bool smoothContinuousFields = args.found("smoothContinuousFields");

    if (lagrangeLimiterTolerance <= 0)
    {
        FatalErrorInFunction
            << "Option -lagrangeLimiterTolerance must be positive, but got "
            << lagrangeLimiterTolerance << exit(FatalError);
    }

    if (kxrcfThreshold <= 0)
    {
        FatalErrorInFunction
            << "Option -kxrcfThreshold must be positive, but got "
            << kxrcfThreshold << exit(FatalError);
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

            exportTimeStep
            (
                runTime,
                dgMesh,
                rho,
                rhoU,
                E,
                thermo(),
                lagrangeLimiterTolerance,
                smoothContinuousFields,
                kxrcfThreshold,
                kxrcfLPRMode
            );
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
