// Assemble the semi-discrete compressible Navier-Stokes-Fourier residual
// cell by cell.
//
// The conservative system is written in divergence/source form as
//
//     dQ/dt + div(Fc(Q) + Fd(Q, gradQ)) = S,
//
// with Q = (rho, rhoU, E).  This readable version builds the fluxes using
// dgExpr algebra on GaussField objects.  Equation-specific assembly lives in
// the three includes below:
//
//     mass     : convective + Franz-Durst mass diffusion
//     momentum : convective + viscous/mass-diffusive stress + gravity
//     energy   : convective + mechanical/thermal fluxes + volume sources
const bool selfDiffusion = thermo->selfDiffusion();

forAll(mesh.cells(), cellI)
{
    // dgBasisField exposes basis values/derivatives for the current cell and
    // is reused by every PDE term assembled on this cell.
    const Foam::dgBasisField& basisField = dgMesh.basisField(cellI);

    GaussField<vector>& UG      = U.gaussFields()[cellI];
    GaussField<scalar>& TG      = T.gaussFields()[cellI];
    GaussField<scalar>& pG      = p.gaussFields()[cellI];
    GaussField<scalar>& EG      = E.gaussFields()[cellI];
    GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];
    GaussField<vector>& SRhoG   = SRho.gaussFields()[cellI];
    GaussField<vector>& SEG     = SE.gaussFields()[cellI];
    GaussField<tensor>& gradUG  = gradU.gaussFields()[cellI];
    GaussField<vector>& gradPG  = gradP.gaussFields()[cellI];
    GaussField<vector>& gradTG  = gradT.gaussFields()[cellI];

    // Dynamic viscosity is cached by dgThermo during update/updateBC from the
    // limited/clamped temperature field. Thermal conductivity is still built
    // on demand from the same stage-consistent T field.
    const GaussField<scalar>& muG = thermo->mu().gaussFields()[cellI];
    GaussField<scalar> kappaG(cellI, &dgMesh);

    thermo->transport().calcKappa(cellI, TG, kappaG);

    #include "stateResiduals/assembleMassEqnResiduals.H"
    #include "stateResiduals/assembleMomentumEqnResiduals.H"
    #include "stateResiduals/assembleEnergyEqnResiduals.H"
}
