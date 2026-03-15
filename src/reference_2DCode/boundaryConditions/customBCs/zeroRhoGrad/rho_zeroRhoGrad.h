#ifndef RHO_ZERORHOGRAD_H
#define RHO_ZERORHOGRAD_H
#include <vector>
namespace zeroRhoGrad
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
    void correctGradRho(std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n, bool isStrongMeth);
}
#endif // RHO_ZERORHOGRAD_H
