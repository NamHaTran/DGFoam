Info<< "Creating conservative fields\n" << endl;

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

