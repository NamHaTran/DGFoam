#include "bcVariables.h"
#include <vector>
#include "../ConstDeclaration.h"

namespace BCVars {
    std::vector<bool> DirichletAppMethUStrong(bcSize,false), NewmannAppMethGradUStrong(bcSize,false),
    DirichletAppMethTStrong(bcSize,false), NewmannAppMethGradTStrong(bcSize,false),
    DirichletAppMethPStrong(bcSize,false), NewmannAppMethGradPStrong(bcSize,false),
    DirichletAppMethGeneralBCStrong(bcSize,false), NewmannAppMethGradGeneralBCStrong(bcSize,false);
}
