// Exchange conservative states across processor boundaries.
// This is typically used right after stage advancement / limiting so the
// processor ghost state matches the authoritative conservative solution.
/*
// Enable this block when you want to inspect plus-side conservative values on
// the tracked processor face before and after synch().
debugProcessorFaceExchange::printFieldPlusValues("before rho.synch()", rho);
debugProcessorFaceExchange::printFieldPlusValues("before rhoU.synch()", rhoU);
debugProcessorFaceExchange::printFieldPlusValues("before E.synch()", E);
*/

rho.synch();
E.synch();
rhoU.synch();

/*
debugProcessorFaceExchange::printFieldPlusValues("after rho.synch()", rho);
debugProcessorFaceExchange::printFieldPlusValues("after rhoU.synch()", rhoU);
debugProcessorFaceExchange::printFieldPlusValues("after E.synch()", E);
*/
