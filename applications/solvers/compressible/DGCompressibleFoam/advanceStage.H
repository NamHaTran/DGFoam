// Apply the explicit stage update to every conservative field once residuals
// are known, then rebuild the primitive/thermo data required by the next stage.
const Switch profileTimings =
    runTime.controlDict().getOrDefault<Switch>("profileTimings", false);

scalar updateStageTime = 0;
scalar limiterTime = 0;
scalar conservativeSynchTime = 0;
scalar primitiveThermoUpdateTime = 0;
scalar thermoSynchTime = 0;
clockTime advanceStageTimer;

if (profileTimings)
{
    advanceStageTimer.resetTime();
}

forAll(mesh.cells(), cellI)
{
    // updateStage() pulls the residual buffer for this cell and applies the
    // active explicit scheme (Euler, RK, SSPRK, ...) to rho/rhoU/E. The time
    // discretization driver selects either global deltaT or local pseudo-deltaT.
    timeDiscretization.updateStage(cellI);
}

if (profileTimings)
{
    updateStageTime = advanceStageTimer.elapsedTime();
    advanceStageTimer.resetTime();
}

// Limit the updated conservative stage state before it is reused by the next
// stage. Limiters act on the conservative DG fields directly, so refresh the
// primitive and thermo fields afterwards.
limiterManager.correct();

if (profileTimings)
{
    limiterTime = advanceStageTimer.elapsedTime();
    advanceStageTimer.resetTime();
}

// The limiter modifies the authoritative conservative stage state, so
// synchronize rho/rhoU/E before rebuilding any primitive variables from it.
#include "synchConservative.H"

if (profileTimings)
{
    conservativeSynchTime = advanceStageTimer.elapsedTime();
    advanceStageTimer.resetTime();
}

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

if (profileTimings)
{
    primitiveThermoUpdateTime = advanceStageTimer.elapsedTime();
    advanceStageTimer.resetTime();
}

// he, a, and mu are thermo-owned fields, so thermo handles their processor
// synchronization internally after the cell-local reconstruction is complete.
thermo->synch();

if (profileTimings)
{
    thermoSynchTime = advanceStageTimer.elapsedTime();
}

// Once the primitive layer has been reconstructed from the synchronized
// conservative state, exchange p/T/U so the next stage sees consistent
// processor-face primitive values.
// #include "synchPrimary.H"

if (profileTimings)
{
    const scalar totalAdvanceStageTime =
        updateStageTime
      + limiterTime
      + conservativeSynchTime
      + primitiveThermoUpdateTime
      + thermoSynchTime;

    Info<< "Profiling advanceStage"
        << ": updateStage=" << updateStageTime
        << " s, limiter=" << limiterTime
        << " s, synchConservative=" << conservativeSynchTime
        << " s, primitiveThermoUpdate=" << primitiveThermoUpdateTime
        << " s, thermoSynch=" << thermoSynchTime
        << " s, total=" << totalAdvanceStageTime
        << " s" << nl;
}
