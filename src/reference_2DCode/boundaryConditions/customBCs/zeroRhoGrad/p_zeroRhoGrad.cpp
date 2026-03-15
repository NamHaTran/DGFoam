#include "p_zeroRhoGrad.h"
#include "VarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/readBCInfor/pressure/readNewmannPresBCValue.h"
#include <fstream>
#include "./boundaryConditions/zeroGradient.h"
#include <vector>

/*Wrap of zeroGradient boundary condition*/
namespace zeroRhoGrad {
    void p_IO(int bcGrp, std::ifstream &FileFlux)
    {
        bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::zeroRhoGrad;

        //Read gradUApp method
        readZeroGradP(FileFlux,bcGrp);
    }

    void correctP(double &pM, double pP)
    {
        //P khong quan trong, de cho co thoi
        pM = zeroGradient_scalar(pP);
    }

    void correctGradP(std::vector<double> &gradM, const std::vector<double> &gradP)
    {
        //P khong quan trong, de cho co thoi
        gradM=gradP;
    }
}
