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

// Primitive values are filled by the thermo reconstruction from rho/rhoU/E.
// Physical boundary ghost states are handled separately by the compressible
// boundary manager, which writes conservative traces and then asks thermo to
// rebuild the matching primitive boundary values.
