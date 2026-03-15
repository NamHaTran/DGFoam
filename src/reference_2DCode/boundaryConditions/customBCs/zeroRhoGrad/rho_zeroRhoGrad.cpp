#include "rho_zeroRhoGrad.h"
#include "./boundaryConditions/zeroGradient.h"
#include <vector>
#include "./boundaryConditions/bcVariables.h"

namespace zeroRhoGrad
{
    /* Cac dieu kien cua rho se overwrite cac dieu kien cua p */
    void correctRho(double &varM, double varP)
    {
        varM = zeroGradient_scalar(varP);
    }

    void correctGradRho(std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n, bool isStrongMeth)
    {

        zeroGradient_correctGrad(gradM,gradP, n, isStrongMeth);
    }
}
