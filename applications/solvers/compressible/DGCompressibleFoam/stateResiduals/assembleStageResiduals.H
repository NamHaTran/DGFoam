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
const scalar selfDiffCoeffScale = thermo->selfDiffCoeffScale();
const dgField<scalar>& muAV =
    artificialDiffusion
  ? artificialDiffusionNonSmooth.updateMuAV()
  : artificialDiffusionNonSmooth.muAV();

forAll(mesh.cells(), cellI)
{
    GaussField<vector>& UG      = U.gaussFields()[cellI];
    GaussField<scalar>& TG      = T.gaussFields()[cellI];
    GaussField<scalar>& pG      = p.gaussFields()[cellI];
    GaussField<scalar>& EG      = E.gaussFields()[cellI];
    GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];
    GaussField<vector>& SRhoG   = SRho.gaussFields()[cellI];
    GaussField<tensor>& SRhoUG  = SRhoU.gaussFields()[cellI];
    GaussField<vector>& SEG     = SE.gaussFields()[cellI];
    GaussField<tensor>& gradUG  = gradU.gaussFields()[cellI];
    GaussField<vector>& gradPG  = gradP.gaussFields()[cellI];
    GaussField<vector>& gradTG  = gradT.gaussFields()[cellI];
    const GaussField<scalar>& muAVG = muAV.gaussFields()[cellI];

    // Transport properties are cached by dgThermo during update/updateBC from
    // the limited/clamped temperature field.
    const GaussField<scalar>& muG = thermo->mu().gaussFields()[cellI];
    const GaussField<scalar>& kappaG =
        thermo->kappa().gaussFields()[cellI];

    tmp<GaussField<vector>> tMassDiffFlux;

    if (selfDiffusion)
    {
        // Franz-Durst mass diffusion is preassembled once per cell and reused across all three PDEs
        const auto massDiffFluxExpr =
           -selfDiffCoeffScale*dg::expr(muG)
           *(
                dg::expr(gradPG)/dg::expr(pG)
              - scalar(0.5)*dg::expr(gradTG)/dg::expr(TG)
            );

        tMassDiffFlux =
            dg::evaluate<vector>(massDiffFluxExpr, cellI, &dgMesh);

        compressibleBC->correctFlux(tMassDiffFlux.ref());
    }

    #include "stateResiduals/assembleMassEqnResiduals.H"
    #include "stateResiduals/assembleMomentumEqnResiduals.H"
    #include "stateResiduals/assembleEnergyEqnResiduals.H"
}
