#include "u_interiorSide.h"
#include "VarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/readBCInfor/velocity/readNewmannVelBCValue.h"
#include <fstream>
#include "./boundaryConditions/zeroGradient.h"
#include <vector>

namespace interiorSide {
    void u_IO(int bcGrp, std::ifstream &FileFlux)
    {
        bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::interiorSide;
        readZeroGradU(FileFlux,bcGrp);
    }
    void correctU(std::vector<double> &varM, std::vector<double> varP)
    {
        zeroGradient_vector(varM, varP);
    }

    void correctGradU(std::vector<double> &gradM, const std::vector<double> &gradP)
    {
        gradM=gradP;
    }
}
