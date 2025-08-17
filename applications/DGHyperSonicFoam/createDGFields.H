Info<< "Reading thermophysical properties\n" << endl;

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Get thermo type and create dgThermo instance
dictionary thermoType(thermoDict.subDict("dgThermo"));
word dgThermoTypeName(thermoType.lookup("type"));

autoPtr<dgThermo> thermo = dgThermo::New(dgThermoTypeName, thermoDict, dgMesh);



//scalar rhoVal = eos->rho(300.0, 101325.0); // T=300K, p=1atm

/*
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);
*/
