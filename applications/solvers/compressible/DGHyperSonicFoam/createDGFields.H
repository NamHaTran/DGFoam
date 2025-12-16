/* Process:
    - Create conservatives fields: rho, rhoU, E
    - Create primary fields: p, T, U
    - Convert conservative to primary (only calculate U)
    - Create thermo: thermo must be created AFTER creating fields because it needs to
    look up the fields to calculate T and p
*/

// Create conservative fields
#include "createConservativeFields.H"

// Create primary fields
#include "createPrimaryFields.H"

// Convert conservative to primary (only calculate U)
#include "conservativeToPrimary.H"

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
// T and p is calculated internally right after themo is constructed

