#ifndef READMIXEDTEMPBCVALUE_H
#define READMIXEDTEMPBCVALUE_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of inletOutlet boundary condition (Mixed type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readInletOutletT(std::ifstream &FileFlux, int bcGrp);
#endif // READMIXEDTEMPBCVALUE_H
