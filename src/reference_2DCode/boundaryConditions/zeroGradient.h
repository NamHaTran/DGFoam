#ifndef ZEROGRADIENT_H
#define ZEROGRADIENT_H
#include <vector>
// zero gradient
void zeroGradient_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP);

double zeroGradient_scalar(double scalarP);

void zeroGradient_correctGrad(std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n, bool isStrongMeth);
#endif // ZEROGRADIENT_H
