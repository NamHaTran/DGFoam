// Advance the explicit time discretization once the stage residuals are ready.
for (label cellI = 0; cellI < mesh.C().size(); cellI++)
{
    timeDiscretization.updateStage(cellI);

    // Rebuild primitive velocity from the updated conservative state.
    U.gaussFields()[cellI] = rhoU.gaussFields()[cellI] / rho.gaussFields()[cellI];

    // Reconstruct all thermodynamic fields from the stage-consistent
    // conservative variables.
    thermo->update(cellI);
}
