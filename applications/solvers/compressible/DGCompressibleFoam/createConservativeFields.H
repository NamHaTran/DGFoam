// Conservative fields are the actual unknowns advanced in time.
// They are read from disk because the solver expects an initial condition
// in terms of rho, rhoU, and total energy E.
Info<< "Create conservative fields\n" << endl;

dgField<scalar> rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    true
);

// rhoU stores the momentum vector, not the primitive velocity.
dgField<vector> rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    true
);

// E denotes the total energy per unit volume used by the conservative
// thermo reconstruction.
dgField<scalar> E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    true
);

// Conservative gradients S = grad(C), with C = (rho, rhoU, E).
// These auxiliary DG fields own DoFs but are initialized to zero and kept out
// of run-time I/O.
dgField<vector> SRho
(
    IOobject
    (
        "SRho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dgMesh,
    true
);

dgField<tensor> SRhoU
(
    IOobject
    (
        "SRhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dgMesh,
    true
);

dgField<vector> SE
(
    IOobject
    (
        "SE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dgMesh,
    true
);
