// Synch data across processor boundaries
scalarProcBC->update(rho);
scalarProcBC->update(p);
scalarProcBC->update(T);
vectorProcBC->update(U);