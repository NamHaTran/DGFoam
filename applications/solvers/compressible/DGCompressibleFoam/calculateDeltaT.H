if (runTime.isAdjustTimeStep())
{
    // Compute a global explicit time step from the current DG state.
    // The estimate is cell-local, then reduced to the minimum over the whole mesh.
    const scalar CFL = runTime.controlDict().getOrDefault<scalar>("maxCo", 0.5);
    const scalar maxDeltaT =
        runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

    if (maxDeltaT <= SMALL)
    {
        FatalIOErrorInFunction(runTime.controlDict())
            << "Entry 'maxDeltaT' must be positive, but got "
            << maxDeltaT << '.'
            << exit(FatalIOError);
    }

    // Explicit DG schemes usually tighten the CFL with the polynomial order.
    // Here we use the common (p+1)^2 scaling.
    const scalar pScale = sqr(scalar(dgMesh.pOrder() + 1));
    scalar dt = GREAT;

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        // h is a geometric characteristic size stored on dgGeomCell.
        const scalar h = dgMesh.cells()[cellI]->cellSize();

        if (h <= SMALL)
        {
            FatalErrorInFunction
                << "Non-positive characteristic cell size for cell " << cellI
                << abort(FatalError);
        }

        const auto lambdaExpr =
            dg::expr(thermo->a().gaussFields()[cellI].cellField())
          + dg::mag(dg::expr(U.gaussFields()[cellI].cellField()));

        // Keep the original CFL semantics: only interior Gauss points
        // contribute to the local time-step estimate.
        const scalar lambdaMax = dg::maxValue(lambdaExpr);

        const scalar dtCell = CFL*h/(pScale*lambdaMax);
        dt = min(dt, dtCell);
    }

    reduce(dt, minOp<scalar>());
    dt = min(dt, maxDeltaT);

    if (dt <= SMALL)
    {
        FatalErrorInFunction
            << "Computed non-positive deltaT = " << dt
            << abort(FatalError);
    }

    runTime.setDeltaT(dt);
}
