#include "rho_reflectRhoGrad.h"
#include "./boundaryConditions/zeroGradient.h"
#include <vector>
#include "./boundaryConditions/bcVariables.h"
#include "DGMath.h"

namespace reflectRhoGrad
{
    /* Cac dieu kien cua rho se overwrite cac dieu kien cua p */
    void correctRho(double &varM, double varP)
    {
        varM = zeroGradient_scalar(varP);
    }

    void correctGradRho(int edgeGrp, std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n)
    {

        double nx(n[0]), ny(n[1]), c(1.0);

        if (
                /* Calculate dot product of grad(rho) and n. If grad(rho) has the same direction with n (dot product > 0), then
                 * mass is compressed at wall and mass diffusion can occur on normal direction pointted to the fluid domain, then
                 * no need to reflect grad(rho) at side Minus (no need to remove normal term of grad(rho)).
                */
                math::vectorDotProduct(gradP,n)>0
                )
        {
            c=0;
        }
        else
        {
            if (BCVars::NewmannAppMethGradGeneralBCStrong[edgeGrp-1])
            {
                /* For reference
                dUXMinus[i] = dUXPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
                dUYMinus[i] = dUYPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;*/

                c=2;
            }
            else
            {
                c=1;
            }
        }

        gradM[0] = gradP[0] - c * (gradP[0] * nx + gradP[1] * ny)*nx;
        gradM[1] = gradP[1] - c * (gradP[0] * nx + gradP[1] * ny)*ny;
    }
}
