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
 * 2. Optionally use positivity-preserving theta fields to clamp selected
 *    conservative variables to DoF0 where requested by
 *    \c smoothByClampingTheta.
 * 3. Apply the internal Lagrange-node positivity limiter to the possibly
 *    clamped fields.
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
    const bool smoothByClampingTheta,
    const scalar clampTheta,
    const volScalarField* theta1FieldPtr,
    const volScalarField* theta2FieldPtr
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

    if (smoothByClampingTheta && theta1FieldPtr && theta2FieldPtr)
    {
        clampConservativeFieldsByTheta
        (
            dgMesh,
            rhoField,
            rhoUField,
            EField,
            clampTheta,
            *theta1FieldPtr,
            *theta2FieldPtr
        );
    }

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

        for (const label pyfrPointI : nodemap)
        {
            const vector xi = pyfrPts[pyfrPointI];
            const vector eta = dgVtkLagrange::xiToEta(type, xi);

            const point x = mapEtaToX(eta, cellVertices, type);

            List<scalar> basis;
            dgVtkLagrange::computeBasisAt(type, eta, pOrder, basis);

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
