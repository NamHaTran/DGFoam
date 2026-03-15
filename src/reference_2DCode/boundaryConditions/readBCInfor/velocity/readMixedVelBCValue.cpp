#include "readMixedVelBCValue.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "../supportReadingBCFuncs.h"
#include "../../fixedValue.h"
#include "../../bcVariables.h"

void readInletOutletU(std::ifstream &FileFlux, int bcGrp)
{
    std::string line, tempStr;

    std::getline(FileFlux, line);
    std::istringstream fixedUStream(line);
    fixedUStream >> tempStr;
    if ((tempStr.compare("inletValue") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'inletValue' in file U/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        fixedUStream >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
    }

    //read U/gradU application method
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong,"U");
    readGradUAppMeth(FileFlux,bcGrp,BCVars::NewmannAppMethGradUStrong,"U");
}
