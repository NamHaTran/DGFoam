// T and p is calculated inside thermo, no need to care
// only calculate U

for (label cellI = 0; cellI < mesh.C().size(); cellI++)
{
    const GaussField<scalar>& rhoG = rho.gaussFields()[cellI]; 
    const GaussField<vector>& rhoUG = rhoU.gaussFields()[cellI];

    GaussField<vector>& UG = U.gaussFields()[cellI];
    
    // Calculate U from rhoU and rho
    UG = rhoUG / rhoG;
}
