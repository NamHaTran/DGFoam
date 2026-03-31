// Exchange primitive states across processor boundaries.
// This is used after rebuilding p, T, and U from the synchronized
// conservative state so boundary reconstruction sees a consistent primitive
// interface state on every rank.
// debugProcessorFaceExchange::printFieldPlusValues("before p.synch()", p);
// debugProcessorFaceExchange::printFieldPlusValues("before T.synch()", T);
// debugProcessorFaceExchange::printFieldPlusValues("before U.synch()", U);

p.synch();
T.synch();
U.synch();

// debugProcessorFaceExchange::printFieldPlusValues("after p.synch()", p);
// debugProcessorFaceExchange::printFieldPlusValues("after T.synch()", T);
// debugProcessorFaceExchange::printFieldPlusValues("after U.synch()", U);
