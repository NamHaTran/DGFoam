/*
Info<< "Reading thermophysical properties\n" << endl;

// Read DOF fields
Info<< "Reading field p\n" << endl;
dofField<scalar> pDof
(
    "p",
    runTime.timeName(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    dgMesh
);

Info<< "Reading field U\n" << endl;
dofField<vector> UDof
(
    "U",
    runTime.timeName(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    dgMesh
);

Info<< "Reading field T\n" << endl;
dofField<scalar> TDof
(
    "T",
    runTime.timeName(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    dgMesh
);

GaussField<scalar> pGauss
(
    &pDof,
    0,
    &dgMesh
);
pGauss.interpolateFromDof();
Info << pGauss << endl;
*/

// Test dgField and GaussField functionalities
dgField<scalar> pField
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    true
);

Info << pField.gaussFields()[0] << endl;

dgField<vector> UField
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    true
);

dgField<scalar> TField
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    true
);

dgField<scalar> TestField
(
    IOobject
    (
        "TestField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    false
);
