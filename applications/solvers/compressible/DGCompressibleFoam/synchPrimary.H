// Exchange primitive states across processor boundaries.
// This is used after rebuilding p, T, and U from the synchronized
// conservative state so boundary reconstruction sees a consistent primitive
// interface state on every rank.
p.synch();
T.synch();
U.synch();
