#ifndef RHO_INTERIORSIDE_H
#define RHO_INTERIORSIDE_H
#include <vector>
namespace interiorSide
{
    /**
     * @brief Function corrects density following interiorSide condition (using zeroGradient_scalar function).
     * @param varM: minus side density
     * @param varP: plus side density
     */
    void correctRho(double &varM, double varP);

    /**
     * @brief Function corrects gradient of density following interiorSide condition.
     *
     * @param gradM: minus side of grad(rho).
     * @param gradP: plus side of grad(rho).
     */
    void correctGradRho(std::vector<double> &gradM, const std::vector<double> &gradP);
}
#endif // RHO_INTERIORSIDE_H
