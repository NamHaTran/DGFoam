#include "readDirichletVelBCValue.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "../supportReadingBCFuncs.h"
#include "../../fixedValue.h"
#include "../../bcVariables.h"

void readNoSlipU(std::ifstream &FileFlux, int bcGrp)
{
    std::string tempStr(""), line;
    bcValues::uBCFixed[bcGrp - 1] = 0.0;
    bcValues::vBCFixed[bcGrp - 1] = 0.0;
    bcValues::wBCFixed[bcGrp - 1] = 0.0;

    //read U/gradU application method
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong,"U");
    //readGradUAppMeth(FileFlux,bcGrp,NewmannAppMethGradUStrong);
}

void readFixedValueU(std::ifstream &FileFlux, int bcGrp)
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
        fixedUStream >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
    }

    //read U/gradU application method
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong,"U");
    //readGradUAppMeth(FileFlux,bcGrp,NewmannAppMethGradUStrong);
}

void readMaxwellSlipU(std::ifstream &FileFlux, int bcGrp)
{
    std::string line, tempStr;
    std::getline(FileFlux, line);
    std::istringstream fixedUStream1(line);
    //Read sigmaU
    fixedUStream1>>tempStr;
    if ((tempStr.compare("sigmaU") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'sigmaU' in file U/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        fixedUStream1>>bcValues::sigmaU;
    }

    std::getline(FileFlux, line);
    std::istringstream fixedUStream2(line);
    //Read UWall
    fixedUStream2 >> tempStr;
    if ((tempStr.compare("uWall") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'uWall' in file U/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        fixedUStream2 >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
    }

    //read U/gradU application method
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong,"U");
    //readGradUAppMeth(FileFlux,bcGrp,NewmannAppMethGradUStrong);
}
