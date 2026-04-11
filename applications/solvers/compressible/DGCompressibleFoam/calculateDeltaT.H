if (runTime.isAdjustTimeStep() || timeDiscretization.localTimeStepping())
{
    // Compute explicit time-step estimates from the current DG state.
    // dgTimeDiscretization decides whether the selected scheme receives the
    // global minimum or capped per-cell pseudo-time increments.
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
    List<scalar> candidateDeltaT(dgMesh.nCells(), runTime.deltaTValue());

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

        if (lambdaMax <= SMALL)
        {
            FatalErrorInFunction
                << "Non-positive maximum wave speed for cell " << cellI
                << ": lambdaMax = " << lambdaMax
                << abort(FatalError);
        }

        const scalar dtCell = min(CFL*h/(pScale*lambdaMax), maxDeltaT);
        candidateDeltaT[cellI] = dtCell;
    }

    timeDiscretization.setDeltaT(candidateDeltaT);
}
