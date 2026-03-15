#include "readSymmetryBC.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "supportReadingBCFuncs.h"
#include "../fixedValue.h"
#include "../bcVariables.h"

void readSymmetryBC(std::ifstream &FileFlux, int bcGrp, std::string fileName)
{
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethGeneralBCStrong,fileName);
    readGradUAppMeth(FileFlux,bcGrp,BCVars::NewmannAppMethGradGeneralBCStrong,fileName);
}
