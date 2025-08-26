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

Info<< "Reading field U\n" << endl;
dgGeneralBoundaryManager<vector> UBC
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading field p\n" << endl;
dgGeneralBoundaryManager<scalar> pBC
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Test boundary conditions
Info<< "\n[TEST] Internal values:\n";
Info<< "U internalField = " << UBC.internalValue() << nl;
Info<< "p internalField = " << pBC.internalValue() << nl;

// Khởi tạo test values
vector testVec(1, 1, 1);
scalar testSca = 5e4;

vector dummyVecGrad(0.123, 0.456, 0.789);
scalar dummyScaGrad = 1.234;

vector plusVec, unusedVecGrad;
scalar plusSca, unusedScaGrad;

Info<< "\n[TEST] Boundary conditions - U\n";
forAll(UBC, patchI)
{
    const dgGeneralBoundaryField<vector>& bc = UBC[patchI];

    // Gọi updateValue
    bc.updateValue
    (
        0,                // gaussID
        vector(1, 0, 0),  // n
        testVec,          // minusValue
        dummyVecGrad,     // minusGrad
        plusVec,          // plusValue
        unusedVecGrad     // plusGrad (unused)
    );

    Info<< "Patch: " << bc.name()
        << ", Type: " << bc.BCtype()
        << ", minus: " << testVec
        << ", plus: " << plusVec << nl;
}

Info<< "\n[TEST] Boundary conditions – p\n";
forAll(pBC, patchI)
{
    const dgGeneralBoundaryField<scalar>& bc = pBC[patchI];

    // Gọi updateValue
    bc.updateValue
    (
        0,
        vector(1, 0, 0),
        testSca,
        dummyScaGrad,
        plusSca,
        unusedScaGrad
    );

    Info<< "Patch: " << bc.name()
        << ", Type: " << bc.BCtype()
        << ", minus: " << testSca
        << ", plus: " << plusSca << nl;
}

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
