#include "readMixedPresBCValue.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "../supportReadingBCFuncs.h"
#include "../../fixedValue.h"
#include "../../bcVariables.h"

void readInletOutletP(std::ifstream &FileFlux, int bcGrp)
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
        fixedUStream >> bcValues::pBCFixed[bcGrp - 1];
    }

    //read U/gradU application method
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethPStrong,"p");
    readGradUAppMeth(FileFlux,bcGrp,BCVars::NewmannAppMethGradPStrong,"p");
}
