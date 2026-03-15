#ifndef READNEWMANNPRESBCVALUE_H
#define READNEWMANNPRESBCVALUE_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of zeroGradient boundary condition (Neumann type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readZeroGradP(std::ifstream &FileFlux, int bcGrp);
#endif // READNEWMANNPRESBCVALUE_H
