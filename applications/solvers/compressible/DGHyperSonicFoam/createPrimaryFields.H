Info<< "Creating primitive fields\n" << endl;

// Primitive fields are created by the solver and updated by thermo.

dgField<scalar> p
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

dgField<scalar> T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    dgMesh,
    false
);

dgField<vector> U
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

// Declare boundary condition managers
IOobject TBCIO
(
    "TBC",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

autoPtr<dgGeneralBoundaryManager<scalar>> TBC
(
    new dgGeneralBoundaryManager<scalar>(TBCIO, dgMesh)
);

IOobject pBCIO
(
    "pBC",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

autoPtr<dgGeneralBoundaryManager<scalar>> pBC
(
    new dgGeneralBoundaryManager<scalar>(pBCIO, dgMesh)
);

IOobject rhoBCIO
(
    "rhoBC",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

autoPtr<dgGeneralBoundaryManager<scalar>> rhoBC
(
    new dgGeneralBoundaryManager<scalar>(rhoBCIO, dgMesh)
);

IOobject UBCIO
(
    "UBC",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

autoPtr<dgGeneralBoundaryManager<vector>> UBC
(
    new dgGeneralBoundaryManager<vector>(UBCIO, dgMesh)
);


// Declare interface boundary condition manager
autoPtr<dgProcessorBoundaryManager<scalar>> scalarProcBC
(
    new dgProcessorBoundaryManager<scalar>(dgMesh)
);

autoPtr<dgProcessorBoundaryManager<vector>> vectorProcBC
(
    new dgProcessorBoundaryManager<vector>(dgMesh)
);
