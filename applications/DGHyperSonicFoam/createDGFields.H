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
