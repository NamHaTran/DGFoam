// dgSchemes collects the DG-specific numerical choices:
// - explicit time scheme,
// - flux scheme per equation term,
// - optional coefficients for each flux/time discretization model.
// This keeps the solver body compact: the tutorial code assembles equations,
// while the exact numerical methods are runtime-selected from dictionaries.
IOdictionary dgSchemesDict
(
    IOobject
    (
        "dgSchemes",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// NonSmooth artificial diffusion owns and registers muAV before the flux
// manager is built, so artificial-diffusion flux solvers can use muAV in their
// penalty coefficients.
dgArtificialDiffusionNonSmooth artificialDiffusionNonSmooth
(
    dgMesh,
    rho,
    U,
    thermo->a()
);

// fluxSolverManager is the indirection layer that picks the right numerical
// flux for each equation term when dgConvectiveTerm asks for a surface flux.
// This lets the PDE assembly code say "mass term", "momentum term", etc.
// without hard-wiring Roe/HLLC/Rusanov logic inside the solver.
dgFluxSolverManager fluxSolverManager
(
    dgSchemesDict,
    dgMesh
);

// limiterManager reads constant/dgOptions and constructs the ordered list of
// DG limiters to be applied after every explicit time stage.
dgLimiterManager limiterManager
(
    dgMesh
);

// timeDiscretization stores residual buffers and applies the chosen explicit
// scheme (Euler, SSPRK, RK4, ...) to all registered conservative fields.
// Conceptually it is the stage driver for the semi-discrete system M dU/dt=R.
dgTimeDiscretization timeDiscretization
(
    dgSchemesDict,
    runTime,
    dgMesh
);

// Only conservative fields are registered because those are the unknowns
// advanced by the explicit update formula. Primitive fields are derived data
// and are reconstructed again after each stage update.
timeDiscretization.registerField(rho);
timeDiscretization.registerField(rhoU);
timeDiscretization.registerField(E);

// Auxiliary mass projection for conservative gradients
//     S = grad(C), C = (rho, rhoU, E).
// The fields are solved algebraically each stage before flux assembly and are
// not time-advanced by dgTimeDiscretization.
dgMassProjection massProjection(dgMesh);

massProjection.registerField(SRho);
massProjection.registerField(SRhoU);
massProjection.registerField(SE);

timeDiscretization.writeInfo(Info);
massProjection.writeInfo(Info);
artificialDiffusionNonSmooth.writeInfo(Info);
limiterManager.listLimiters(Info);
