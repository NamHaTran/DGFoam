#ifndef RHO_REFLECTRHOGRAD_H
#define RHO_REFLECTRHOGRAD_H
#include <vector>
namespace reflectRhoGrad
{
    /**
     * @brief Function corrects density following reflectRhoGrad condition (using zeroGradient_scalar function).
     * @param pM: minus side density
     * @param pP: plus side density
     */
    void correctRho(double &varM, double varP);

    /**
     * @brief Function corrects gradient of density following reflectRhoGrad condition.
     *
     * Remove normal gradient component.
     *
     * @param gradM: minus side grad(rho)
     * @param gradP: plus side grad(rho)
     */
    void correctGradRho(int edgeGrp, std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n);
}
#endif // RHO_REFLECTRHOGRAD_H
