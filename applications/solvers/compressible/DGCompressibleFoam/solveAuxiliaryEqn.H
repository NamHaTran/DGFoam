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
//    auxEqnResidual assembles the weak gradient equation
//
//        M S = - int_Omega C grad(phi) dV
//              + int_dOmega phi C* n dS
//
//    where C is rho, rhoU, or E.  Internal and processor faces use the
//    LDG-style central state C* = 0.5*(C- + C+).  Physical boundary faces use
//    the middle state already supplied by updateBC.H through
//    faceGaussField::midValues().
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

dgGradProjectionTerm<scalar> rhoGradProjectionTerm
(
    "rhoGradProjectionTerm",
    dgMesh
);

dgGradProjectionTerm<scalar> eGradProjectionTerm
(
    "eGradProjectionTerm",
    dgMesh
);

dgGradProjectionTerm<vector> rhoUGradProjectionTerm
(
    "rhoUGradProjectionTerm",
    dgMesh
);

forAll(mesh.cells(), cellI)
{
    // Step 1a: gather the stage-consistent conservative Gauss fields for
    // this cell.  prepareStageFields.H has already synchronized primitive
    // values and applied conservative boundary conditions, so face minus,
    // plus, and boundary mid states are ready for the auxiliary projection.
    GaussField<scalar>& rhoG  = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG = rhoU.gaussFields()[cellI];
    GaussField<scalar>& EG    = E.gaussFields()[cellI];

    // Step 1b: register one conservative state field per projection term and
    // calculate the central face flux stored on the term's face field.
    rhoGradProjectionTerm.add(rhoG);
    rhoGradProjectionTerm.calculateFlux(cellI);

    eGradProjectionTerm.add(EG);
    eGradProjectionTerm.calculateFlux(cellI);

    rhoUGradProjectionTerm.add(rhoUG);
    rhoUGradProjectionTerm.calculateFlux(cellI);

    auxEqnResidual<scalar> scalarAuxResidual(dgMesh, cellI);
    const label rhoAuxTermI =
        scalarAuxResidual.addTerm(rhoGradProjectionTerm);
    const label eAuxTermI =
        scalarAuxResidual.addTerm(eGradProjectionTerm);

    auxEqnResidual<vector> vectorAuxResidual(dgMesh, cellI);
    const label rhoUAuxTermI =
        vectorAuxResidual.addTerm(rhoUGradProjectionTerm);

    // Step 1c: assemble the modal residual for M S = R. Physical boundary
    // C* values are read from the term flux; internal and processor C* values
    // were computed as 0.5*(C- + C+) by calculateFlux().
    scalarAuxResidual.calcResidual();
    vectorAuxResidual.calcResidual();

    // Step 1d: store the three residual vectors in the auxiliary projection
    // driver.  The actual division by the mass-matrix diagonal is done after
    // all cells have contributed.
    massProjection.addResidual
    (
        SRho.name(),
        cellI,
        scalarAuxResidual.R(rhoAuxTermI)
    );

    massProjection.addResidual
    (
        SRhoU.name(),
        cellI,
        vectorAuxResidual.R(rhoUAuxTermI)
    );

    massProjection.addResidual
    (
        SE.name(),
        cellI,
        scalarAuxResidual.R(eAuxTermI)
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
SRho.interpolateFromDof();
SRhoU.interpolateFromDof();
SE.interpolateFromDof();

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
