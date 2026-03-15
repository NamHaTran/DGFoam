#include "supportReadingBCFuncs.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>

void readUAppMeth(std::ifstream &FileFlux, int bcGrp, std::vector<bool> &UMethod, std::string fileName)
{
    std::string tempStr(""), line;
    std::getline(FileFlux, line);
    std::istringstream stream(line);

    //Read UAppMethod (application method)
    stream >> tempStr;
    if ((tempStr.compare("UAppMethod") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'UAppMethod' in file "+fileName+"/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        stream >> tempStr;
        if ((tempStr.compare("weak") == 0))
        {
            UMethod[bcGrp - 1] = false;

        }
        else if ((tempStr.compare("strong") == 0))
        {
            UMethod[bcGrp - 1] = true;
        }
        else
        {
            message::writeLog(systemVar::pwd, systemVar::caseName, "Unknow method " +tempStr+ " of UAppMethod in file "+fileName+"/group " + std::to_string(bcGrp) + ".\n");
        }
    }
}

void readGradUAppMeth(std::ifstream &FileFlux, int bcGrp, std::vector<bool> &gradUMethod, std::string fileName)
{
    std::string tempStr(""), line;
    std::getline(FileFlux, line);
    std::istringstream stream(line);

    stream >> tempStr;
    if ((tempStr.compare("gradUAppMethod") != 0))
    {
        message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'gradUAppMethod' in file "+fileName+"/group " + std::to_string(bcGrp) + ".\n");
    }
    else
    {
        stream >> tempStr;
        if ((tempStr.compare("weak") == 0))
        {
            gradUMethod[bcGrp - 1] = false;

        }
        else if ((tempStr.compare("strong") == 0))
        {
            gradUMethod[bcGrp - 1] = true;
        }
        else
        {
            message::writeLog(systemVar::pwd, systemVar::caseName, "Unknow method " +tempStr+ " of gradUAppMethod in file "+fileName+"/group " + std::to_string(bcGrp) + ".\n");
        }
    }
}
