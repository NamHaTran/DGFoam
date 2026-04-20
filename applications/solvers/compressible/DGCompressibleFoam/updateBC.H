// For the current owner cell:
// 1. impose ghost conservative boundary data on rho/rhoU/E,
// 2. reconstruct primitive and thermo traces from that updated conservative
//    boundary state.
//
// This makes the compressible boundary manager the single source of truth for
// physical boundary conditions while still feeding U/p/T to the flux routines.
GaussField<vector>& UG    = U.gaussFields()[cellI];
GaussField<scalar>& EG    = E.gaussFields()[cellI];
GaussField<scalar>& rhoG  = rho.gaussFields()[cellI];
GaussField<vector>& rhoUG = rhoU.gaussFields()[cellI];

compressibleBC->update(rhoG, rhoUG, EG, UG);

// Rebuild p/T/he/a/mu on the just-updated ghost conservative state so the
// boundary face traces used by DG fluxes are fully consistent.
thermo->updateBC(cellI);
