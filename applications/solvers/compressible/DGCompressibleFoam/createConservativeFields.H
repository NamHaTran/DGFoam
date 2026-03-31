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
