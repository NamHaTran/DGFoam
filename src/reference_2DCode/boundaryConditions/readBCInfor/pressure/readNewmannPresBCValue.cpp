#include "readNewmannPresBCValue.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "../supportReadingBCFuncs.h"
#include "../../fixedValue.h"
#include "../../bcVariables.h"

void readZeroGradP(std::ifstream &FileFlux, int bcGrp)
{
    //readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong);
    readGradUAppMeth(FileFlux,bcGrp,BCVars::NewmannAppMethGradPStrong,"p");
    //Set gia tri tam cho dk bien
    bcValues::pBCFixed[bcGrp - 1]=iniValues::pIni;
}
