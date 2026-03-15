#include "fixedValue.h"
#include <vector>
#include "zeroGradient.h"

void weakFixedValue_vector(std::vector<double> &vectorM, const std::vector<double> &constVector)
{
    vectorM=constVector;
}

void strongFixedValue_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP, const std::vector<double> &constVector)
{
    //zeroGradient_vector(vectorM, constVector);
    //Voi vector, constVector = 0.5*(vectorP+vectorM)
    for (int i=0; i<static_cast<int>(vectorP.size()); i++)
    {
        vectorM[i]=2*constVector[i]-vectorP[i];
    }
}

void fixedValue_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP, const std::vector<double> &constVector, bool isStrongMeth)
{
    if (isStrongMeth)
        strongFixedValue_vector(vectorM,vectorP,constVector);
    else
        weakFixedValue_vector(vectorM,constVector);
}

double weakFixedValue_scalar(double constScalar)
{
    return constScalar;
}

double strongFixedValue_scalar(double scalarP, double constScalar)
{
    return (2*constScalar-scalarP);
}

double fixedValue_scalar(double scalarP, double constScalar, bool isStrongMeth)
{
    if (isStrongMeth)
        return strongFixedValue_scalar(scalarP,constScalar);
    else
        return weakFixedValue_scalar(constScalar);
}
