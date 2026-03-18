// Exchange primitive states across processor boundaries.
// This happens before boundary reconstruction so remote neighboring cells
// contribute up-to-date interface data on parallel runs.
// The synchronized objects are dgField values, so this communication updates
// the DG face data needed to reconstruct correct interface states.
scalarProcBC->update(p);
scalarProcBC->update(T);
vectorProcBC->update(U);
