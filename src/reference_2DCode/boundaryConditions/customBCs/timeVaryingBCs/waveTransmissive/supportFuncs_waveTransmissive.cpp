#include "supportFuncs_waveTransmissive.h"
#include "VarDeclaration.h"
#include <fstream>
#include <sstream>
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/readBCInfor/supportReadingBCFuncs.h"
#include <vector>

namespace waveTransmissive
{
    bool *includeSoundSpeed_p = new bool[bcSize];
    bool *includeSoundSpeed_T = new bool[bcSize];
    bool *includeSoundSpeed_u = new bool[bcSize];

    void readCondition(int bcGrp, std::ifstream &FileFlux, std::string file)
    {
        if ((file.compare("U") == 0))
            bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::waveTransmissive;
        else if ((file.compare("p") == 0))
            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::waveTransmissive;
        else if ((file.compare("T") == 0))
            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::waveTransmissive;


        std::string line, tempStr, tempStr2;
        std::getline(FileFlux, line);
        std::istringstream fixedUStream1(line);
        //Read soundSpeed
        fixedUStream1>>tempStr;
        if ((tempStr.compare("soundSpeed") != 0))
        {
            message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'soundSpeed' in file "+file+"/group " + std::to_string(bcGrp) + ".\n");
        }
        else
        {
            fixedUStream1>>tempStr2;

            if ((file.compare("U") == 0))
            {
                if ((tempStr2.compare("yes") == 0)||(tempStr2.compare("true") == 0))
                    waveTransmissive::includeSoundSpeed_u[bcGrp - 1]=true;
                else
                    waveTransmissive::includeSoundSpeed_u[bcGrp - 1]=false;

                //Set gia tri ban dau cho dk bien
                bcValues::uBCFixed[bcGrp - 1]=0;
                bcValues::uBCFixed[bcGrp - 1]=0;
            }
            else if ((file.compare("p") == 0))
            {
                if ((tempStr2.compare("yes") == 0)||(tempStr2.compare("true") == 0))
                    waveTransmissive::includeSoundSpeed_p[bcGrp - 1]=true;
                else
                    waveTransmissive::includeSoundSpeed_p[bcGrp - 1]=false;

                //Set gia tri ban dau cho dk bien
                bcValues::pBCFixed[bcGrp - 1]=iniValues::pIni;
            }
            else if ((file.compare("T") == 0))
            {
                if ((tempStr2.compare("yes") == 0)||(tempStr2.compare("true") == 0))
                    waveTransmissive::includeSoundSpeed_T[bcGrp - 1]=true;
                else
                    waveTransmissive::includeSoundSpeed_T[bcGrp - 1]=false;

                //Set gia tri ban dau cho dk bien
                bcValues::TBCFixed[bcGrp - 1]=iniValues::TIni;
            }
        }

        //read U/gradU application method
        readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong,file);
    }
}
