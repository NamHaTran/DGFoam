#ifndef RHO_H
#define RHO_H
#include <vector>
namespace zeroRhoGradUncorectP
{
    void correctRho(double &varM, double varP);

    void correctGradRho(std::vector<double> &gradM, const std::vector<double> &gradP);
}
#endif // RHO_H
