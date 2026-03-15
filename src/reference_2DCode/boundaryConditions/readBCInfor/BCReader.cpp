#include "./BCReader.h"
#include "dynamicVarDeclaration.h"
#include <iostream>
#include <string>

//Boundary condition libs
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/fixedValue.h"

#include "./boundaryConditions/readBCInfor/velocity/readDirichletVelBCValue.h"
#include "./boundaryConditions/readBCInfor/velocity/readNewmannVelBCValue.h"
#include "./boundaryConditions/readBCInfor/velocity/readMixedVelBCValue.h"

#include "./boundaryConditions/readBCInfor/pressure/readNewmannPresBCValue.h"
#include "./boundaryConditions/readBCInfor/pressure/readDirichletPresBCValue.h"
#include "./boundaryConditions/readBCInfor/pressure/readMixedPresBCValue.h"

#include "./boundaryConditions/readBCInfor/temperature/readNewmannTempBCValue.h"
#include "./boundaryConditions/readBCInfor/temperature/readDirichletTempBCValue.h"
#include "./boundaryConditions/readBCInfor/temperature/readMixedTempBCValue.h"

#include "./boundaryConditions/readBCInfor/readSymmetryBC.h"


//Custom Boundary Conditions -------------------------------------------------
//Key word de tim trong file: Custom_Boundary_Conditions
//Dieu kien bien zeroRhoGradUncorectP
#include "./boundaryConditions/customBCs/zeroRhoGradUncorectP/p_zeroRhoGradUncorrectP.h"

//Dieu kien bien reflect
#include "./boundaryConditions/customBCs/reflectRhoGrad/p_reflectRhoGrad.h"

//Dieu kien bien interiorSide
#include "./boundaryConditions/customBCs/interiorSide/p_interiorSide.h"
#include "./boundaryConditions/customBCs/interiorSide/T_interiorSide.h"
#include "./boundaryConditions/customBCs/interiorSide/u_interiorSide.h"

//Dieu kien bien zeroRhoGrad
#include "./boundaryConditions/customBCs/zeroRhoGrad/p_zeroRhoGrad.h"

//Dieu kien bien Maxwell
#include "./boundaryConditions/customBCs/nonEquilibriumBCs/MaxwellSlip/u_MaxwellSlip.h"
#include "./boundaryConditions/customBCs/nonEquilibriumBCs/SmoluchowskyTJump/T_SmoluchowskyTJump.h"

//Dieu kien bien waveTransmissive
#include "./boundaryConditions/customBCs/timeVaryingBCs/waveTransmissive/supportFuncs_waveTransmissive.h"
//Custom Boundary Conditions -------------------------------------------------

namespace readVectorBC {
    void u(std::string mode)
    {
        /* Velocity boundary conditions:
            Id = 1:---------------------------------------------------------------------------------------------------------
                fixedValue
                value u v w
                Description: ham nay set dieu kien fixedValue cho patch.

            Id = 2:---------------------------------------------------------------------------------------------------------
                noSlip

            Id = 3:---------------------------------------------------------------------------------------------------------
                movingWall
                value u v w
                Description: set velocity cho wall (typically: cavity case)

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                inletValue u v w
                Description: Neu subsonic fixedValue neu inFlow va zeroGradient neu outFlow. Neu supersonic deu la zeroGradient.
                Nen dung cho outlet.

            Id = 5:---------------------------------------------------------------------------------------------------------
                slip
                sigmaU value
                uWall u v w
                Description: apply dieu kien slip cua Maxwell

            Id = 6:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry

            Id = 10:--------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song
        */

        std::string fileName("U.txt"), tempStr(""), Loc;
        if (mode.compare("p")==0)
        {
            Loc = (systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + std::to_string(systemVar::currentProc) + "/0");
        }
        else {
            Loc = (systemVar::wD + "/CASES/" + systemVar::caseName + "/0");
        }

        if (systemVar::currentProc==0)
        {
            std::cout << "	Reading " << fileName << "\n";
        }

        std::string FileLoc(Loc + "/" + fileName);
        std::ifstream FileFlux(FileLoc.c_str());
        int bcGrp(0), bcIndex(0);
        std::vector<bool> check_bc(meshVar::nBc,false);

        if (FileFlux)
        {
            std::string line, keyWord;
            while (std::getline(FileFlux, line))
            {
                std::istringstream line2str(line);
                std::vector<std::string> ptr;
                //Split <line2str> stringstream into array of words
                while ((line2str >> keyWord))
                {
                    ptr.push_back(keyWord);
                }

                int numWrd = static_cast<int>(ptr.size());
                if (numWrd >= 2 && ptr[0].compare("//"))
                {
                    std::string str1(ptr[0]);

                    if (bcIndex<meshVar::nBc)
                    {
                        for (int i=0; i<meshVar::nBc; i++)
                        {
                            if (!check_bc[i])
                            {
                                if (str1.compare("initialValue") == 0)  //initial data
                                {
                                    std::istringstream str_u(ptr[1]);
                                    std::istringstream str_v(ptr[2]);
                                    std::istringstream str_w(ptr[3]);
                                    str_u >> iniValues::uIni;
                                    str_v >> iniValues::vIni;
                                    str_w >> iniValues::wIni;
                                    goto label;
                                }
                                else if (str1.compare("Type") == 0)  //Bc Type
                                {
                                    std::string str0(ptr[1]);

                                    if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
                                    {
                                        if ((str0.compare("noSlip") == 0))  //Type noSlip
                                        {
                                            bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::noSlip;
                                            readNoSlipU(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("slip") == 0))  //Type slip
                                        {
                                            //Use Maxwell-Smoluchovsky boundary condition
                                            bcValues::slipBCFlag=true;
                                            bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::slipWall;
                                            readMaxwellSlipU(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("movingWall") == 0))  //Type movingWall
                                        {
                                            bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::movingWall;
                                            readFixedValueU(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//Custom_Boundary_Conditions----------------------------------------
                                        //19/08/2021 MaxwellSlip
                                        else if ((str0.compare("MaxwellSlip") == 0))
                                        {
                                            MaxwellSlip::u_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//------------------------------------------------------------------
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "wall"));
                                        }
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
                                    {
                                        if ((str0.compare("fixedValue") == 0))  //Type fixedValue
                                        {
                                            bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::fixedValue;
                                            readFixedValueU(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("inletOutlet") == 0))  //Type inletOutlet
                                        {
                                            bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::inletOutlet;
                                            readInletOutletU(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }

                                        // 11/09/2020: add them zeroGradient cho velocityBC, Id = 6
                                        else if ((str0.compare("zeroGradient") == 0))  //Type zeroGradient
                                        {
                                            bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::zeroGrad;
                                            readZeroGradU(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//Custom_Boundary_Conditions----------------------------------------
                                        //22/07/2021 interiorSide
                                        else if ((str0.compare("interiorSide") == 0))
                                        {
                                            interiorSide::u_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }


                                        //03/06/2022 waveTransmissive
                                        else if ((str0.compare("waveTransmissive") == 0))
                                        {
                                            bcValues::timeVaryingBCAvailable=true;
                                            waveTransmissive::readCondition(bcGrp,FileFlux,"U");

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//------------------------------------------------------------------

                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "patch"));
                                        }
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)  //SYMMETRY
                                    {
                                        if ((str0.compare("symmetry") == 0))  //Type symmetry
                                        {
                                            bcValues::UBcType[bcGrp - 1] = BCVars::generalBCId::symmetry;
                                            readSymmetryBC(FileFlux,bcGrp,"U");

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "symmetry"));
                                        }
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 4)  //MATCHED
                                    {
                                        if ((str0.compare("matched") == 0))  //Type matched
                                        {
                                            bcValues::UBcType[bcGrp - 1] = BCVars::generalBCId::matched;
                                            bcValues::uBCFixed[bcGrp - 1] = 0.0;
                                            bcValues::vBCFixed[bcGrp - 1] = 0.0;
                                            bcValues::wBCFixed[bcGrp - 1] = 0.0;
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "matched"));
                                        }
                                    }
                                }
                                else if ((str1.compare("Group") == 0))  //Group
                                {
                                    std::istringstream str_bcGrp(ptr[1]);
                                    str_bcGrp >> bcGrp;
                                }
                            }
                        }
                    }
                }
                label:
                if ((bcIndex>=meshVar::nBc))
                {
                    break;
                }
            }
        }
        else
        {
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
        }
    }
}

namespace readScalarBC {
    void p(std::string mode)
    {
        /* Pressure boundary conditions:
            Id = 1:---------------------------------------------------------------------------------------------------------
                fixedValue
                value p

            Id = 2:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid. Argument pP de nguyen, de phong truong hop
                muon lay zeroGradient la gia tri tai diem Gauss dang tinh cua cell phia plus.

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                inletValue p
                Description:

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry

            Id = 10:---------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song

            Id = 11:---------------------------------------------------------------------------------------------------------
                interpFrmDensity
                Description: Dung cho wall. Gia tri cua P duoc tinh tu rhoP va Twall, van dung zeroGradient cho P khi tinh gradient.
                correctRho = true khi dung kieu BC nay.
        */

        std::string fileName = "p.txt";
        std::string tempStr("");
        std::string Loc;
        int bcIndex(0);
        std::vector<bool> check_bc(meshVar::nBc,false);
        if (mode.compare("p")==0)
        {
            Loc=systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + std::to_string(systemVar::currentProc) + "/0";
        }
        else {
            Loc=systemVar::wD + "/CASES/" + systemVar::caseName + "/0";
        }

        if (systemVar::currentProc==0)
        {
            std::cout << "	Reading " << fileName << "\n";
        }
        std::string FileLoc(Loc + "/" + fileName);
        std::ifstream FileFlux(FileLoc.c_str());

        int bcGrp(0);

        if (FileFlux)
        {
            std::string line, keyWord;
            while (std::getline(FileFlux, line))
            {
                std::istringstream line2str(line);
                std::vector<std::string> ptr;
                //Split <line2str> stringstream into array of words
                while ((line2str >> keyWord))
                {
                    ptr.push_back(keyWord);
                }

                int numWrd = static_cast<int>(ptr.size());
                if (numWrd >= 2 && ptr[0].compare("//"))
                {
                    std::string str1(ptr[0]);
                    if (bcIndex<meshVar::nBc)
                    {
                        for (int i=0; i<meshVar::nBc; i++)
                        {
                            if (!check_bc[i])
                            {
                                if (str1.compare("initialValue") == 0)  //initial data
                                {
                                    std::istringstream str_val(ptr[1]);
                                    str_val >> iniValues::pIni;
                                    goto label;
                                }
                                else if (str1.compare("Type") == 0)  //Bc Type
                                {
                                    std::string str0(ptr[1]);

                                    if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
                                    {
                                        if ((str0.compare("zeroGradient") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::zeroGrad;
                                            readZeroGradP(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("fixedValue") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::fixedValue;
                                            readFixedValueP(FileFlux,bcGrp);
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("inletOutlet") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::inletOutlet;
                                            readInletOutletP(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("interpFrmDensity") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::interpFrmDensity;
                                            readZeroGradP(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }

//------------------------------------------//Custom_Boundary_Conditions----------------------------------------
                                        //zeroGradRhoUncorrectP
                                        else if ((str0.compare("zeroGradRhoUncorrectP") == 0))
                                        {
                                            zeroRhoGradUncorectP::p_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//------------------------------------------------------------------

                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "wall/patch"));
                                        }
                                        check_bc[i]=true;
                                        bcIndex++;
                                        goto label;
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
                                    {
                                        if ((str0.compare("zeroGradient") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::zeroGrad;
                                            readZeroGradP(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("fixedValue") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::fixedValue;
                                            readFixedValueP(FileFlux,bcGrp);
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("inletOutlet") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::inletOutlet;
                                            readInletOutletP(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("interpFrmDensity") == 0))
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::pressureBCId::interpFrmDensity;
                                            readZeroGradP(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//Custom_Boundary_Conditions----------------------------------------
                                        //22/07/2021 interiorSide
                                        else if ((str0.compare("interiorSide") == 0))
                                        {
                                            interiorSide::p_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }

                                        //01/07/2021 reflectGradRho
                                        else if ((str0.compare("reflectGradRho") == 0))
                                        {
                                            reflectRhoGrad::p_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }

                                        //23/09/2021 zeroRhoGrad (for patch)
                                        else if ((str0.compare("zeroRhoGrad") == 0))
                                        {
                                            zeroRhoGrad::p_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }

                                        //03/06/2022 waveTransmissive (for patch)
                                        else if ((str0.compare("waveTransmissive") == 0))
                                        {
                                            bcValues::timeVaryingBCAvailable=true;
                                            waveTransmissive::readCondition(bcGrp,FileFlux,"p");

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//------------------------------------------------------------------
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)
                                    {
                                        if ((str0.compare("symmetry") == 0))  //Type symmetry
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::generalBCId::symmetry;

                                            readSymmetryBC(FileFlux,bcGrp,"p");

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "symmetry"));
                                        }
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 4)
                                    {
                                        if ((str0.compare("matched") == 0))  //Type matched
                                        {
                                            bcValues::pBcType[bcGrp - 1] = BCVars::generalBCId::matched;
                                            bcValues::pBCFixed[bcGrp - 1] = iniValues::pIni;
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "matched"));
                                        }
                                    }
                                }
                                else if ((str1.compare("Group") == 0))  //Group
                                {
                                    std::istringstream str_bcGrp(ptr[1]);
                                    str_bcGrp >> bcGrp;
                                }
                            }
                        }
                    }
                }
                label:
                if ((bcIndex>=meshVar::nBc))
                {
                    break;
                }
            }
        }
        else
        {
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
        }
    }

    void T(std::string mode)
    {
        /* Temperature boundary conditions:
            Id = 2:---------------------------------------------------------------------------------------------------------
                fixedValue
                value T
                ==> grad surface = grad surface plus

            Id = 3:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid. Argument TP de nguyen, de phong truong hop
                muon lay zeroGradient la gia tri tai diem Gauss dang tinh cua cell phia plus.
                ==> giai zero normal grad

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                inletValue T
                Description:
                ==> giai zero normal grad neu outflow, grad surface = grad mean neu inflow

            Id = 6:---------------------------------------------------------------------------------------------------------
                temperatureJump
                sigmaT value
                TWall T
                Description: apply dieu kien temperature cua Smoluchowsky
                ==> giai zero normal grad

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry
                ==> reflect vector grad T

            Id = 10:---------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song
        */

        std::string fileName = "T.txt";
        std::string tempStr("");
        std::string Loc;
        if (mode.compare("p")==0)
        {
            Loc=systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + std::to_string(systemVar::currentProc) + "/0";
        }
        else {
            Loc=systemVar::wD + "/CASES/" + systemVar::caseName + "/0";
        }

        if (systemVar::currentProc==0)
        {
            std::cout << "	Reading " << fileName << "\n";
        }
        std::string FileLoc(Loc + "/" + fileName);
        std::ifstream FileFlux(FileLoc.c_str());
        int bcIndex(0);
        std::vector<bool> check_bc(meshVar::nBc,false);

        int bcGrp(0);

        if (FileFlux)
        {
            std::string line, keyWord;
            while (std::getline(FileFlux, line))
            {
                std::istringstream line2str(line);
                std::vector<std::string> ptr;
                //Split <line2str> stringstream into array of words
                while ((line2str >> keyWord))
                {
                    ptr.push_back(keyWord);
                }

                int numWrd = static_cast<int>(ptr.size());
                if (numWrd >= 2 && ptr[0].compare("//"))
                {
                    std::string str1(ptr[0]);
                    if (bcIndex<meshVar::nBc)
                    {
                        for (int i=0; i<meshVar::nBc; i++)
                        {
                            if (!check_bc[i])
                            {
                                if (str1.compare("initialValue") == 0)  //initial data
                                {
                                    std::istringstream str_val(ptr[1]);
                                    str_val >> iniValues::TIni;
                                    goto label;
                                }
                                else if (str1.compare("Type") == 0)  //Bc Type
                                {
                                    std::string str0(ptr[1]);

                                    if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
                                    {
                                        if ((str0.compare("zeroGradient") == 0))
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::zeroGrad;
                                            readZeroGradT(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("fixedValue") == 0))
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::fixedValue;
                                            readFixedValueT(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("inletOutlet") == 0))
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::inletOutlet;
                                            readInletOutletT(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("temperatureJump") == 0))  //Type temperatureJump
                                        {
                                            bcValues::temperatureJump=true;
                                            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::temperatureJump;
                                            readTemperatureJump(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//Custom_Boundary_Conditions----------------------------------------
                                        //22/07/2021 interiorSide
                                        else if ((str0.compare("interiorSide") == 0))
                                        {
                                            interiorSide::T_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        //19/08/2021 interiorSide
                                        else if ((str0.compare("SmoluchowskyTJump") == 0))
                                        {
                                            SmoluchowskyTJump::T_IO(bcGrp,FileFlux);

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//------------------------------------------------------------------
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "wall/patch"));
                                        }
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
                                    {
                                        if ((str0.compare("zeroGradient") == 0))
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::zeroGrad;
                                            readZeroGradT(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("fixedValue") == 0))
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::fixedValue;
                                            readFixedValueT(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
                                        else if ((str0.compare("inletOutlet") == 0))
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::inletOutlet;
                                            readInletOutletT(FileFlux,bcGrp);

                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//Custom_Boundary_Conditions----------------------------------------
                                        //03/06/2022 waveTransmissive (for patch)
                                        else if ((str0.compare("waveTransmissive") == 0))
                                        {
                                            bcValues::timeVaryingBCAvailable=true;
                                            waveTransmissive::readCondition(bcGrp,FileFlux,"T");

                                            //Luon co 3 dong code nay
                                            check_bc[i]=true;
                                            bcIndex++;
                                            goto label;
                                        }
//------------------------------------------//------------------------------------------------------------------
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)
                                    {
                                        if ((str0.compare("symmetry") == 0))  //Type symmetry
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::generalBCId::symmetry;
                                            readSymmetryBC(FileFlux,bcGrp,"T");
                                        }
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "symmetry"));
                                        }
                                        check_bc[i]=true;
                                        bcIndex++;
                                        goto label;
                                    }
                                    else if (meshVar::BoundaryType[bcGrp - 1][1] == 4)
                                    {
                                        if ((str0.compare("matched") == 0))  //Type matched
                                        {
                                            bcValues::TBcType[bcGrp - 1] = BCVars::generalBCId::matched;
                                            bcValues::TBCFixed[bcGrp - 1] = iniValues::TIni;
                                        }
                                        else
                                        {
                                            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "matched"));
                                        }
                                        check_bc[i]=true;
                                        bcIndex++;
                                        goto label;
                                    }
                                }
                                else if ((str1.compare("Group") == 0))  //Group
                                {
                                    std::istringstream str_bcGrp(ptr[1]);
                                    str_bcGrp >> bcGrp;
                                }
                            }
                        }
                    }
                }
                label:
                if ((bcIndex>=meshVar::nBc))
                {
                    break;
                }
            }
        }
        else
        {
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
        }

        for (int i = 0; i < meshVar::nBc; i++)
        {
            if (bcValues::TBcType[i]!=6 && (bcValues::UBcType[i]==5))
            {
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::SlipBcCompatibleError(i+1));
            }
        }
    }
}
