#include "readNewmannVelBCValue.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "../supportReadingBCFuncs.h"
#include "../../fixedValue.h"
#include "../../bcVariables.h"

void readZeroGradU(std::ifstream &FileFlux, int bcGrp)
{
    //readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong);
    readGradUAppMeth(FileFlux,bcGrp,BCVars::NewmannAppMethGradUStrong,"U");
    //Set gia tri tam cho dk bien
    bcValues::uBCFixed[bcGrp - 1]=0;
    bcValues::vBCFixed[bcGrp - 1]=0;
}
