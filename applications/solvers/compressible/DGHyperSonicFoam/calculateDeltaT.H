if (runTime.isAdjustTimeStep())
{
    // Compute a global explicit time step from the current DG state.
    // The estimate is cell-local, then reduced to the minimum over the whole mesh.
    const scalar CFL = runTime.controlDict().getOrDefault<scalar>("maxCo", 0.5);

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

        // The time-step constraint is evaluated at interior Gauss points so the
        // estimate sees the local high-order DG state, not just one cell average.
        // cellGaussField behaves like a small per-cell field object, and the math
        // layer overloads operators on it pointwise.
        const cellGaussField<vector>& UCell = U.gaussFields()[cellI].cellField();
        const cellGaussField<scalar>& aCell = thermo->a().gaussFields()[cellI].cellField();

        // lambdaMax approximates the stiffest local wave/transport speed.
        // Build the pointwise sum in smaller tmp-backed steps to avoid
        // over-referencing the same temporary cellGaussField object.
        tmp<cellGaussField<scalar>> tLambdaCell = aCell + mag(UCell);

        const scalar lambdaMax = maxValue(tLambdaCell);

        const scalar dtCell = CFL*h/(pScale*lambdaMax);
        dt = min(dt, dtCell);
    }

    reduce(dt, minOp<scalar>());

    if (dt <= SMALL)
    {
        FatalErrorInFunction
            << "Computed non-positive deltaT = " << dt
            << abort(FatalError);
    }

    runTime.setDeltaT(dt);
}
