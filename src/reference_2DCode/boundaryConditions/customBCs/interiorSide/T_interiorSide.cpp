#include "T_interiorSide.h"
#include "VarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/readBCInfor/temperature/readNewmannTempBCValue.h"
#include <fstream>
#include "./boundaryConditions/zeroGradient.h"
#include <vector>

namespace interiorSide {
    void T_IO(int bcGrp, std::ifstream &FileFlux)
    {
        bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::interiorSide;
        readZeroGradT(FileFlux,bcGrp);
    }
    void correctT(double &varM, double varP)
    {
        varM = zeroGradient_scalar(varP);
    }

    void correctGradT(std::vector<double> &gradM, const std::vector<double> &gradP)
    {
        gradM=gradP;
    }
}
