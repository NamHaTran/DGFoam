#ifndef READMIXEDVELBCVALUE_H
#define READMIXEDVELBCVALUE_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of inletOutlet boundary condition (Mixed type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readInletOutletU(std::ifstream &FileFlux, int bcGrp);
#endif // READMIXEDVELBCVALUE_H
