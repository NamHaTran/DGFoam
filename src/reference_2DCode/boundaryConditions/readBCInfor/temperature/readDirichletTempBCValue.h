#ifndef READDIRICHLETTEMPBCVALUE_H
#define READDIRICHLETTEMPBCVALUE_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of fixedValue boundary condition (Dirichlet type).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readFixedValueT(std::ifstream &FileFlux, int bcGrp);
/**
 * @brief Function read inputs of Temperature Jump boundary condition (Dirichlet type) (OLD).
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readTemperatureJump(std::ifstream &FileFlux, int bcGrp);
#endif // READDIRICHLETTEMPBCVALUE_H
