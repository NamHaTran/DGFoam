// Prepare all data needed before residual assembly of the current explicit stage.
// This is the "make stage input state consistent" step.
// A good mental model is:
// - dgTimeDiscretization owns the stage state,
// - the solver makes all derived Gauss fields consistent with that state,
// - only then can dgGeneralPDETerm assemble residuals safely.

// Boundary reconstruction uses primitive interface values, so start each
// stage from a synchronized p/T/U processor state.
#include "synchPrimary.H"

forAll(mesh.cells(), cellI)
{
    // Residuals are recomputed from scratch at every stage.
    timeDiscretization.clearResiduals(cellI);

    // Apply primitive BCs and rebuild the corresponding ghost conservative
    // states so each cell face has a complete left/right state for fluxes.
    // In DG terms, this fills the face Gauss points seen by the owner cell.
    #include "updateBC.H"
}

// Flux solvers cache face states; clear them after the stage input state has
// changed and before new residuals are assembled.
#include "resetState.H"
