// Solve the LDG auxiliary-gradient block for the current explicit stage.
//
// The viscous Navier-Stokes-Fourier residual needs grad(U) and grad(T), but
// the time integrator advances conservative variables.  This include therefore
// does the gradient update in two layers:
//
// 1. Project conservative gradients with the local DG mass matrix:
//
//        SRho  = grad(rho)
//        SRhoU = grad(rhoU)
//        SE    = grad(E)
//
//    dgGradProjectionTerm assembles the weak gradient equation
//
//        M S = - int_Omega C grad(phi) dV
//              + int_dOmega phi C* n dS
//
//    where C is rho, rhoU, or E.  Internal faces use the LDG-style central
//    state C* = 0.5*(C- + C+).  Physical boundary faces use the middle state
//    already supplied by updateBC.H through faceGaussField::midValues().
//
// 2. Rebuild primitive/thermal gradients from those conservative gradients.
//    In particular, gradU is not solved as an independent DoF field here; it
//    is reconstructed pointwise from SRho and SRhoU by the quotient rule:
//
//        grad(U) = (grad(rhoU) - U grad(rho))/rho.
//
//    thermo->updateGradient() then applies the energy/EOS chain rules to get
//    grad(T) and grad(p) from SRho, SE, and gradU.

massProjection.clearResiduals();

forAll(mesh.cells(), cellI)
{
    // Step 1a: gather the stage-consistent conservative Gauss fields for
    // this cell.  prepareStageFields.H has already synchronized primitive
    // values and applied conservative boundary conditions, so face minus,
    // plus, and boundary mid states are ready for the auxiliary projection.
    const Foam::dgBasisField& basisField = dgMesh.basisField(cellI);

    const GaussField<scalar>& rhoG  = rho.gaussFields()[cellI];
    const GaussField<vector>& rhoUG = rhoU.gaussFields()[cellI];
    const GaussField<scalar>& EG    = E.gaussFields()[cellI];

    // Step 1b: create one local LDG projection term per conservative state.
    // Each term owns a local copy of its GaussField because updateFlux()
    // writes the central internal-face state into midValues().
    dgGradProjectionTerm<scalar> auxGradRhoTerm
    (
        "gradRhoFlux",
        basisField,
        rhoG
    );

    dgGradProjectionTerm<vector> auxGradRhoUTerm
    (
        "gradRhoUFlux",
        basisField,
        rhoUG
    );

    dgGradProjectionTerm<scalar> auxGradETerm
    (
        "gradEFlux",
        basisField,
        EG
    );

    // Step 1c: fill C* on internal faces, then assemble the modal residual
    // for M S = R.  Boundary C* values are intentionally left as provided by
    // the boundary-condition path.
    auxGradRhoTerm.updateFlux();
    auxGradRhoUTerm.updateFlux();
    auxGradETerm.updateFlux();

    auxGradRhoTerm.calcResidual();
    auxGradRhoUTerm.calcResidual();
    auxGradETerm.calcResidual();

    // Step 1d: store the three residual vectors in the auxiliary projection
    // driver.  The actual division by the mass-matrix diagonal is done after
    // all cells have contributed.
    massProjection.addResidual
    (
        SRho.name(),
        cellI,
        auxGradRhoTerm.R()
    );

    massProjection.addResidual
    (
        SRhoU.name(),
        cellI,
        auxGradRhoUTerm.R()
    );

    massProjection.addResidual
    (
        SE.name(),
        cellI,
        auxGradETerm.R()
    );
}

forAll(mesh.cells(), cellI)
{
    // Step 1e: solve the cell-local mass projections.  dgMassProjection writes
    // the modal DoFs of SRho, SRhoU, and SE, then interpolates the owner-cell
    // Gauss values for this cell.
    massProjection.solve(cellI);
}

// Refresh all Gauss values once more after every cell DoF has been solved.
// This is needed because face plus/minus values depend on neighbouring cell
// DoFs, which may have been updated after an earlier cell was projected.
forAll(mesh.cells(), cellI)
{
    SRho.gaussFields()[cellI].interpolateFromDof();
    SRhoU.gaussFields()[cellI].interpolateFromDof();
    SE.gaussFields()[cellI].interpolateFromDof();
}

SRho.synch();
SRhoU.synch();
SE.synch();

forAll(mesh.cells(), cellI)
{
    // Step 2a: reconstruct the primitive velocity gradient from the projected
    // conservative gradients.  This keeps the auxiliary solve conservative,
    // while the viscous stress below can still use grad(U) directly.
    const GaussField<scalar>& rhoG = rho.gaussFields()[cellI];
    const GaussField<vector>& UG = U.gaussFields()[cellI];
    const GaussField<vector>& SRhoG = SRho.gaussFields()[cellI];
    const GaussField<tensor>& SRhoUG = SRhoU.gaussFields()[cellI];

    GaussField<tensor>& gradUG = gradU.gaussFields()[cellI];

    dg::assign
    (
        gradUG,
        (
            dg::expr(SRhoUG)
          - dg::expr(UG)*dg::expr(SRhoG)
        )/dg::expr(rhoG)
    );
}

gradU.synch();

forAll(mesh.cells(), cellI)
{
    // Step 2b: let thermo finish the chain rules:
    // - grad(he) from rho, E, U, SRho, SE, gradU
    // - grad(T) from the selected energy/EOS model
    // - grad(p) from the selected EOS
    thermo->updateGradient(cellI);
}

thermo->synchGradient();

massProjection.writeResiduals(Info);

massProjection.clearResiduals();

// Gradient-aware solvers may cache face-local state. The auxiliary projection
// has changed gradient fields, so make the following residual assembly fresh.
#include "resetState.H"
