// Apply the explicit stage update to every conservative field once residuals
// are known, then rebuild the primitive/thermo data required by the next stage.
forAll(mesh.cells(), cellI)
{
    // updateStage() pulls the residual buffer for this cell and applies the
    // active explicit scheme (Euler, RK, SSPRK, ...) to rho/rhoU/E. The time
    // discretization driver selects either global deltaT or local pseudo-deltaT.
    timeDiscretization.updateStage(cellI);
}

// Limit the updated conservative stage state before it is reused by the next
// stage. Limiters act on the conservative DG fields directly, so refresh the
// primitive and thermo fields afterwards.
limiterManager.correct();

// The limiter modifies the authoritative conservative stage state, so
// synchronize rho/rhoU/E before rebuilding any primitive variables from it.
#include "synchConservative.H"

forAll(mesh.cells(), cellI)
{
    // Rebuild the kinematic primitive field outside thermo so the thermo layer
    // remains focused on thermal variables.
    const GaussField<scalar>& rhoG = rho.gaussFields()[cellI];
    const GaussField<vector>& rhoUG = rhoU.gaussFields()[cellI];
    GaussField<vector>& UG = U.gaussFields()[cellI];

    dg::assign
    (
        UG.cellField(),
        dg::expr(rhoUG.cellField())/dg::expr(rhoG.cellField())
    );

    dg::assign
    (
        UG.faceField(),
        dg::expr(rhoUG.faceField())/dg::expr(rhoG.faceField())
    );

    // Reconstruct p, T, he, a, gamma, etc. from the stage-consistent
    // conservative state. This keeps the primitive and thermo layers in sync
    // before the next stage begins.
    thermo->update(cellI);
}

// he and a are thermo-owned fields, so thermo handles their processor
// synchronization internally after the cell-local reconstruction is complete.
// debugProcessorFaceExchange::printFieldPlusValues
// (
//     "before thermo->synch() he",
//     thermo->he()
// );
// debugProcessorFaceExchange::printFieldPlusValues
// (
//     "before thermo->synch() a",
//     thermo->a()
// );

thermo->synch();

// debugProcessorFaceExchange::printFieldPlusValues
// (
//     "after thermo->synch() he",
//     thermo->he()
// );
// debugProcessorFaceExchange::printFieldPlusValues
// (
//     "after thermo->synch() a",
//     thermo->a()
// );

// Once the primitive layer has been reconstructed from the synchronized
// conservative state, exchange p/T/U so the next stage sees consistent
// processor-face primitive values.
// #include "synchPrimary.H"
