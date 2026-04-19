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

    // The transport properties are evaluated through the selected thermo
    // transport model, so Sutherland/power-law/etc. stay runtime-selected.
    // muG and kappaG are full GaussField objects because the diffusive flux
    // needs both cell values and face minus/plus values.
    GaussField<scalar> muG(cellI, &dgMesh);
    GaussField<scalar> kappaG(cellI, &dgMesh);

    thermo->transport().calcMu(cellI, TG, muG);
    thermo->transport().calcKappa(cellI, TG, kappaG);

    #include "stateResiduals/assembleMassEqnResiduals.H"
    #include "stateResiduals/assembleMomentumEqnResiduals.H"
    #include "stateResiduals/assembleEnergyEqnResiduals.H"
}
