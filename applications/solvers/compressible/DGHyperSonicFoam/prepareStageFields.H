// Prepare primitive, conservative, and thermo fields for the current stage.
for (label cellI = 0; cellI < mesh.C().size(); cellI++)
{
    timeDiscretization.clearResiduals(cellI);

    GaussField<vector>& UG      = U.gaussFields()[cellI];
    GaussField<scalar>& TG      = T.gaussFields()[cellI];
    GaussField<scalar>& EG      = E.gaussFields()[cellI];
    GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];
    GaussField<vector>& rhoUG   = rhoU.gaussFields()[cellI];

    // Impose boundary data on the primitive fields used to reconstruct the
    // conservative state of the current RK stage.
    TBC->updateValue(TG);
    rhoBC->updateValue(rhoG);
    UBC->updateValue(UG);

    rhoUG = rhoG*UG;

    GaussField<scalar> eG(cellI, &dgMesh);
    thermo->thermo().calcInternalE(cellI, TG, eG);

    tmp<GaussField<scalar>> tK = 0.5*magSqr(UG);
    EG = rhoG*(eG + tK());

    tK.clear();
}
