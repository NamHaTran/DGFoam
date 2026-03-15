#ifndef READNEWMANNVELBCVALUE_H
#define READNEWMANNVELBCVALUE_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of zeroGradient boundary condition (Neumann type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readZeroGradU(std::ifstream &FileFlux, int bcGrp);
#endif // READNEWMANNVELBCVALUE_H
