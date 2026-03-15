#ifndef READMIXEDPRESBCVALUE_H
#define READMIXEDPRESBCVALUE_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of inletOutlet boundary condition (Mixed type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readInletOutletP(std::ifstream &FileFlux, int bcGrp);
#endif // READMIXEDPRESBCVALUE_H
