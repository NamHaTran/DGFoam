#include "readDirichletTempBCValue.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "../supportReadingBCFuncs.h"
#include "../../fixedValue.h"
#include "../../bcVariables.h"

void readFixedValueT(std::ifstream &FileFlux, int bcGrp)
{
    std::string line, tempStr;
    std::getline(FileFlux, line);
    std::istringstream fixedUStream(line);
    fixedUStream >> tempStr;
    if ((tempStr.compare("value") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'value' in file U/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        fixedUStream >> bcValues::TBCFixed[bcGrp - 1];
    }

    //read U/gradU application method
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethTStrong,"T");
    //readGradUAppMeth(FileFlux,bcGrp,NewmannAppMethGradUStrong);
}

void readTemperatureJump(std::ifstream &FileFlux, int bcGrp)
{
    std::string line, tempStr;
    std::getline(FileFlux, line);
    std::istringstream Stream1(line);
    //Read sigmaT
    Stream1>> tempStr;
    if ((tempStr.compare("sigmaT") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'sigmaT' in file T/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        Stream1>> bcValues::sigmaT;
    }

    std::getline(FileFlux, line);
    std::istringstream Stream2(line);
    //Read TWall
    Stream2 >> tempStr;
    if ((tempStr.compare("TWall") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'TWall' in file T/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        Stream2 >> bcValues::TBCFixed[bcGrp - 1];
    }

    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethTStrong,"T");
}
