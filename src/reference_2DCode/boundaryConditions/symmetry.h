#ifndef SYMMETRY_H
#define SYMMETRY_H
#include "vector"
#include "ConstDeclaration.h"
//void symmetry_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP,  const std::vector<double> &n);
//double symmetry_scalar(double scalarP);

/* Control flag of strong/weak type
 * Strong/weak type defines how to treat normal terms at minus side
 * - strong: normTerm- = -normTerm+: this type gives reflection effect faster but is unstable in high velocity/gradient case (can test them)
 * - weak: normTerm- = 0: this type gives reflection effect slower but is more stable in high velocity/gradient case
*/

void symmetryBC_auxEqs(std::vector <std::vector<double>> &Fluxes, int element, int edge, int edgeGrp, int nG);
void symmetryBC_NSFEqs(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG);
#endif // SYMMETRY_H
