#ifndef READDIRICHLETVELBCVALUE_H
#define READDIRICHLETVELBCVALUE_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of fixedValue boundary condition (Dirichlet type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readNoSlipU(std::ifstream &FileFlux, int bcGrp);
void readFixedValueU(std::ifstream &FileFlux, int bcGrp);
void readMaxwellSlipU(std::ifstream &FileFlux, int bcGrp);
#endif // READDIRICHLETVELBCVALUE_H
