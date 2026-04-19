// Exchange conservative states across processor boundaries.
// This is typically used right after stage advancement / limiting so the
// processor ghost state matches the authoritative conservative solution.
rho.synch();
E.synch();
rhoU.synch();
