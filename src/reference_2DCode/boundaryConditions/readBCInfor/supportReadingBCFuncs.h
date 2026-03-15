#ifndef SUPPORTREADINGBCFUNCS_H
#define SUPPORTREADINGBCFUNCS_H
#include <sstream>
#include <vector>
/**
 * @brief Function reads application method of field p, T, U.
 *
 * There are 2 method: strong and weak.
 *
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 * @param UMethod: flag of application method (returned value).
 */
void readUAppMeth(std::ifstream &FileFlux, int bcGrp, std::vector<bool> &UMethod, std::string fileName);
/**
 * @brief Function reads application method of field grad(p, T, U).
 *
 * There are 2 method: strong and weak.
 *
 * @param FileFlux: file flux to read.
 * @param bcGrp: group ID of boundary.
 * @param UMethod: flag of application method (returned value).
 */
void readGradUAppMeth(std::ifstream &FileFlux, int bcGrp, std::vector<bool> &gradUMethod, std::string fileName);

#endif // SUPPORTREADINGBCFUNCS_H
