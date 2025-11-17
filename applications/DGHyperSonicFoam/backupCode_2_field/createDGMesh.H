Info<< "Creating DG Geometric Mesh with polymonial order" << pOrder << "\n" << endl;

dgGeomMesh dgMesh 
(
    mesh,
    pOrder
);

/*
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


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

// Read DOF of conservative fields
// Helper lambda: check file existence
auto fieldExists = [&](const word& fieldName)
{
    fileName fPath = runTime.path() / runTime.timeName() / fieldName;
    return isFile(fPath);
};

// Loop DOFs
for (label d = 0; d < nDof; ++d)
{
    word rhoName  = word("rho")  + Foam::name(d);
    word rhoUName = word("rhoU") + Foam::name(d);
    word rhoEName = word("rhoE") + Foam::name(d);

    if (fieldExists(rhoName) && fieldExists(rhoUName) && fieldExists(rhoEName))
    {
        // === Field exists → read it ===
        rhoDOFs[d] = volScalarField
        (
            IOobject
            (
                rhoName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        rhoUDOFS[d] = volVectorField
        (
            IOobject
            (
                rhoUName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        rhoEDOFs[d] = volScalarField
        (
            IOobject
            (
                rhoEName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        // === Field missing → initialize zero ===
        rhoDOFs[d] = volScalarField
        (
            IOobject
            (
                rhoName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimDensity, 0.0)
        );

        rhoUDOFS[d] = volVectorField
        (
            IOobject
            (
                rhoUName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimDensity*dimVelocity, vector::zero)
        );

        rhoEDOFs[d] = volScalarField
        (
            IOobject
            (
                rhoEName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimDensity*dimEnergy, 0.0)
        );
    }
}

// Initialize the first DOF with initial conditions
if (runTime.timeName() == "0")
{
    rhoDOFs[0]  = rhoInit;
    rhoUDOFS[0] = rhoDOFs[0] * U;
    rhoEDOFs[0] = rhoDOFs[0] * (Cv*T + 0.5*magSqr(U));
}
*/