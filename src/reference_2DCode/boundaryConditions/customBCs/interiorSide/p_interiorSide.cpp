#include "p_interiorSide.h"
#include "VarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/readBCInfor/pressure/readNewmannPresBCValue.h"
#include <fstream>
#include "./boundaryConditions/zeroGradient.h"
#include <vector>

namespace interiorSide {
    void p_IO(int bcGrp, std::ifstream &FileFlux)
    {
        bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::interiorSide;
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
