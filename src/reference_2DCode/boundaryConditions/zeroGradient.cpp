#include "zeroGradient.h"
#include <vector>
#include <tuple>
#include <DGAuxUltilitiesLib.h>

// zero gradient
void zeroGradient_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP)
{
    vectorM=vectorP;
}

double zeroGradient_scalar(double scalarP)
{
    return (scalarP);
}

void zeroGradient_correctGrad(std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n, bool isStrongMeth)
{
    /* Nguyen tac: normal grad minus = 0
     *             tangential grad minus = tangential grad plus
    */

    if (isStrongMeth)
    {
        gradM[0] = gradP[0] - 2 * (gradP[0] * n[0] + gradP[1] * n[1])*n[0];
        gradM[1] = gradP[1] - 2 * (gradP[0] * n[0] + gradP[1] * n[1])*n[1];
    }
    else
    {
        double tangGrad(-gradP[0]*n[1]+gradP[1]*n[0]);
        gradM[0]=-tangGrad*n[1];
        gradM[1]=tangGrad*n[0];
    }

}

