// Primitive fields are reconstructed from the conservative variables and
// are used to evaluate fluxes and boundary conditions in a more natural form.
Info<< "Create primitive fields\n" << endl;

// The solver owns p, T, and U, but their values are refreshed by the thermo
// model and by stage-update logic after rho/rhoU/E change. These are still
// dgField objects, so they carry high-order Gauss-point data, not only one
// cell-average value per control volume.

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

// Boundary managers read user-specified DG boundary data from constant/
// and know how to inject those values into face Gauss points. The result is
// a boundary state defined at quadrature points, which is exactly what the
// numerical flux routines need.
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

// Initialize primary fields by initial values read from BC settings in constant
TBC->initializeField(T);
pBC->initializeField(p);
UBC->initializeField(U);

// Processor boundary managers synchronize primitive fields across partition
// interfaces before each residual assembly on parallel runs. They play the
// same role as processor patches in OpenFOAM, but at the DG Gauss-point level.
autoPtr<dgProcessorBoundaryManager<scalar>> scalarProcBC
(
    new dgProcessorBoundaryManager<scalar>(dgMesh)
);

autoPtr<dgProcessorBoundaryManager<vector>> vectorProcBC
(
    new dgProcessorBoundaryManager<vector>(dgMesh)
);
