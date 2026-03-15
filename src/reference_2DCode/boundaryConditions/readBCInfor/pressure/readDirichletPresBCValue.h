#ifndef READDIRICHLETPRESBCVALUE_H
#define READDIRICHLETPRESBCVALUE_H
#include <sstream>
#include <vector>

/**
 * @brief Function read inputs of fixedValue boundary condition (Dirichlet type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readFixedValueP(std::ifstream &FileFlux, int bcGrp);
#endif // READDIRICHLETPRESBCVALUE_H
