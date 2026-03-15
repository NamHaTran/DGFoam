#ifndef FIXEDVALUE_H
#define FIXEDVALUE_H
#include <vector>
#include "ConstDeclaration.h"

void weakFixedValue_vector(std::vector<double> &vectorM, const std::vector<double> &constVector);

void strongFixedValue_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP, const std::vector<double> &constVector);

void fixedValue_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP, const std::vector<double> &constVector, bool isStrongMeth);

double weakFixedValue_scalar(double constScalar);

double strongFixedValue_scalar(double scalarP, double constScalar);

double fixedValue_scalar(double scalarP, double constScalar, bool isStrongMeth);

#endif // FIXEDVALUE_H
