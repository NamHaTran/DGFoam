// Synch only the stage input primitive fields across processor boundaries.
scalarProcBC->update(rho);
scalarProcBC->update(T);
vectorProcBC->update(U);
