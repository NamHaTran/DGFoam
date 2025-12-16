
// Read solver settings for DGHyperSonicFoam
IOdictionary fluxSolverDict
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

dgFluxSolverManager fluxSolverManager
(
    fluxSolverDict,
    dgMesh
);
