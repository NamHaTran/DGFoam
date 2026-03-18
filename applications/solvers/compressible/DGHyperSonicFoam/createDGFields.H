/*
Tutorial note:
- The conservative fields must exist first because they are the solved
  variables and the thermo model looks them up by name from the registry.
- The primitive fields must also exist before thermo construction because
  thermo writes reconstructed p/T/U into those DG fields.
- Once the thermo object is constructed, it can update all primitive and
  transport quantities cell by cell from rho, rhoU, and E.
- A dgField can be thought of as a whole-mesh DG variable; each cell entry is
  a GaussField, and each GaussField contains both cellGaussField and
  faceGaussField data. That layering is what makes local DG algebra compact.
*/

// Create conservative fields
#include "createConservativeFields.H"

// Create primary fields
#include "createPrimaryFields.H"

// Thermo reconstruction is separated from the solver loop so that the solver
// can ask for p, T, a, mu, gamma, etc. without duplicating thermodynamic
// formulas in every include file. In practice, thermo becomes the bridge
// between conservative unknowns and the primitive quantities used in fluxes.
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

// The constructor immediately looks up the already-created DG fields and
// initializes thermo-derived quantities on all cells.
autoPtr<dgThermoConservative> thermo = dgThermoConservative::New(dgThermoTypeName, thermoDict, dgMesh);
