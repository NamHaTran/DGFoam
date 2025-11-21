Info<< "Creating conservatives fields\n" << endl;

dgField<scalar> rhoField
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

dgField<vector> rhoUField
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

dgField<scalar> EField
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

Info<< "Creating primary fields\n" << endl;

dgField<scalar> pField
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    false
);

dgField<vector> UField
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    false
);


// Declare thermo
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

autoPtr<dgThermoConservative> thermo = dgThermoConservative::New(dgThermoTypeName, thermoDict, dgMesh);

