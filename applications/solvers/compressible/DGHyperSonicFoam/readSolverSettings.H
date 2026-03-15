
// Read solver settings for DGHyperSonicFoam
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

dgFluxSolverManager fluxSolverManager
(
    dgSchemesDict,
    dgMesh
);

dgTimeDiscretization timeDiscretization
(
    dgSchemesDict,
    runTime,
    dgMesh
);

timeDiscretization.registerField(rho);
timeDiscretization.registerField(rhoU);
timeDiscretization.registerField(E);

timeDiscretization.writeInfo(Info);
