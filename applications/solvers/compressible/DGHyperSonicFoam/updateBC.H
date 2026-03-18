// For the current owner cell:
// 1. impose primitive boundary data on the face Gauss points,
// 2. reconstruct ghost conservative variables from those primitive values,
// 3. update the remaining thermo quantities from that already-imposed
//    boundary state.
//
// The ordering matters: rhoBasedConservative::updateBC() assumes the primary
// ghost values (U, T, p) are already authoritative and must not recompute
// them internally.
// The key data structure here is boundaryGaussField: a lightweight container
// holding only boundary-face quadrature values. It lets us do boundary algebra
// with the same scalar/vector syntax used for interior Gauss fields.
GaussField<vector>& UG      = U.gaussFields()[cellI];
GaussField<scalar>& pG      = p.gaussFields()[cellI];
GaussField<scalar>& TG      = T.gaussFields()[cellI];
GaussField<scalar>& EG      = E.gaussFields()[cellI];
GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];

TBC->updateValue(TG);
pBC->updateValue(pG);
UBC->updateValue(UG);

// Extract the just-imposed primitive ghost states from the face field.
// extractBCGhostState() flattens the boundary quadrature values into a
// boundaryGaussField so they can be transformed independently of geometry.
boundaryGaussField<scalar> pb = pG.faceField().extractBCGhostState();
boundaryGaussField<scalar> Tb = TG.faceField().extractBCGhostState();
boundaryGaussField<vector> Ub = UG.faceField().extractBCGhostState();

// Convert primitive BC data into conservative BC data so the DG flux routines
// can work with the same variable set everywhere.
// Notice how the code reads almost like the formulas themselves:
//     rhoU_b = rho_b * U_b
//     E_b    = rho_b * (e_b + 0.5*|U_b|^2)
boundaryGaussField<scalar> rhob(pb.size());
thermo->eos().calcRhoFromPT(pb, Tb, rhob);
rhoG.faceField().assignBCGhostState(rhob);

boundaryGaussField<vector> rhoUb = rhob*Ub;
rhoUG.faceField().assignBCGhostState(rhoUb);

boundaryGaussField<scalar> eb(Tb.size());
thermo->thermo().calcInternalE(Tb, eb);

boundaryGaussField<scalar> Eb = rhob*(eb + 0.5*magSqr(Ub));
EG.faceField().assignBCGhostState(Eb);

// Fill the derived thermo boundary fields (a, mu, gamma, etc.) while keeping
// the primitive ghost values imposed above as the authoritative BC state.
thermo->updateBC(cellI);
