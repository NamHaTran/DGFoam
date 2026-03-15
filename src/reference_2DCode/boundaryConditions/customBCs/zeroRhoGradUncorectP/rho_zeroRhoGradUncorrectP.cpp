#include "rho_zeroRhoGradUncorrectP.h"
#include "./boundaryConditions/zeroGradient.h"
#include <vector>
namespace zeroRhoGradUncorectP
{
    /* Cac dieu kien cua rho se overwrite cac dieu kien cua p */
    void correctRho(double &varM, double varP)
    {
        varM = zeroGradient_scalar(varP);
    }

    void correctGradRho(std::vector<double> &gradM, const std::vector<double> &gradP)
    {
        gradM=gradP;
    }
}
