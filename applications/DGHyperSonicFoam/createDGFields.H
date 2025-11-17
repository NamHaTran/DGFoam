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

