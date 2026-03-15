#ifndef MATCHED_H
#define MATCHED_H
#include "vector"
void matchedBC_auxEqs(std::vector <std::vector<double>> &Fluxes, int edge, int nG);

void matchedBC_NSFEqs(std::vector<double> &Fluxes, int element, int edge, int nG);
#endif // MATCHED_H
