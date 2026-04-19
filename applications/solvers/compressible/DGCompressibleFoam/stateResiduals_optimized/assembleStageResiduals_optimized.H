// Optimized semi-discrete compressible residual assembly.
//
// This version mirrors stateResiduals/assembleStageResiduals.H but builds
// algebraic flux/source fields with explicit loops over Gauss points.
forAll(mesh.cells(), cellI)
{
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

    GaussField<scalar> muG(cellI, &dgMesh);
    GaussField<scalar> kappaG(cellI, &dgMesh);

    thermo->transport().calcMu(cellI, TG, muG);
    thermo->transport().calcKappa(cellI, TG, kappaG);

    #include "stateResiduals_optimized/assembleMassEqnResiduals_optimized.H"
    #include "stateResiduals_optimized/assembleMomentumEqnResiduals_optimized.H"
    #include "stateResiduals_optimized/assembleEnergyEqnResiduals_optimized.H"
}
