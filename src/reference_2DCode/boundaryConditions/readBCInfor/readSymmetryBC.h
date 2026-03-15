#ifndef READSYMMETRYBC_H
#define READSYMMETRYBC_H
#include <sstream>
#include <vector>
/**
 * @brief Function read inputs of symmetry boundary condition.
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 */
void readSymmetryBC(std::ifstream &FileFlux, int bcGrp, std::string fileName);
#endif // READSYMMETRYBC_H
