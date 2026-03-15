#include "nonEqmBCs_GenFuncs.h"
#include "DGAuxUltilitiesLib.h"
#include "VarDeclaration.h"
#include "DGIOLib.h"
#include "dynamicVarDeclaration.h"
#include <iostream>
#include <sstream>

#include "../../bcVariables.h"

//Time-varying boundary conditions
#include "MaxwellSlip/u_MaxwellSlip.h"
#include "SmoluchowskyTJump/T_SmoluchowskyTJump.h"

namespace nonEquilibriumBCs {
    void updateBCs()
    {
        /*
         * Ham update gia tri tren cac field cua surfaceBCFields, chi dung cho cac BC bien thien theo thoi gian.
         * Hien tai, ham su dung cho temperatureJump va slip conditions
        */
        if (auxUlti::checkNonEqmBCAvailable() && flowProperties::viscous)
        {
            if (systemVar::currentProc==0)
            {
                std::cout<<"Updating nonequilibrium BCs.\n";
            }

            int globleEdge(0);
            for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
            {
                globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
                int edgeGrp(auxUlti::getGrpOfEdge(globleEdge));
                int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]);

                for (int nG=0; nG<=mathVar::nGauss1D; nG++)
                {
                    //Update Temperature
                    if (TType == BCVars::temperatureBCId::SmoluchowskyTJump)
                    {
                        //SmoluchowskyTJump::calcTJump_DGTypeExplicit(globleEdge,edgeGrp,nG);
                        SmoluchowskyTJump::calcTJump_FDMTypeSemiImplicit(globleEdge,edgeGrp,nG);
                        //SmoluchowskyTJump::calcTJump_FDMTypeImplicit(globleEdge,edgeGrp,nG);
                    }

                    //Update Velocity
                    if (UType == BCVars::velocityBCId::MaxwellSlip)
                    {
                        //MaxwellSlip::calcUSlip_DGTypeExplicit(globleEdge,edgeGrp,nG);
                        MaxwellSlip::calcUSlip_FDMTypeImplicit(globleEdge,edgeGrp,nG);
                    }
                }
            }
        }
    }

    void resizeSurfaceFields()
    {
        meshVar::GaussPtsOnBCEdge_x = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,0.0);
        meshVar::GaussPtsOnBCEdge_y = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,0.0);
        meshVar::GaussPtsOnBCEdge_unitVector_x = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,0.0);
        meshVar::GaussPtsOnBCEdge_unitVector_y = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,0.0);
        meshVar::distanceFromGaussPtsToCentroid = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,0.0);

        /*
        if (auxUlti::checkNonEqmBCAvailable())
        {
            SurfaceBCFields::TBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,iniValues::TIni);
            SurfaceBCFields::uBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,iniValues::uIni);
            SurfaceBCFields::vBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,iniValues::vIni);
            SurfaceBCFields::pBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss1D+1,iniValues::pIni);
        }*/
    }
}

namespace nonEquilibriumBCs_IO {
    void readSurfaceValues(std::string Loc)
    {
        if (auxUlti::checkNonEqmBCAvailable())
        {
            std::string fileLoc(""), fileName("");

            //Read file TSurface & USurface
            fileName = "TSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::readFromFile(fileLoc,fileName,false,SurfaceBCFields::TBc);
            fileName = "uSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::readFromFile(fileLoc,fileName,false,SurfaceBCFields::uBc);
            fileName = "vSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::readFromFile(fileLoc,fileName,false,SurfaceBCFields::vBc);
        }
    }

    void writeSurfaceValues(std::string Loc)
    {
        if (auxUlti::checkNonEqmBCAvailable())
        {
            std::string fileLoc(""), fileName("");

            //Write file TSurface & USurface
            fileName = "TSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::writeToFile(fileLoc,fileName,SurfaceBCFields::TBc);
            fileName = "uSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::writeToFile(fileLoc,fileName,SurfaceBCFields::uBc);
            fileName = "vSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::writeToFile(fileLoc,fileName,SurfaceBCFields::vBc);
        }
    }

    void writeToFile(std::string loc, std::string name, double **array)
    {
        std::ofstream Flux(loc.c_str());
        int numRow(meshVar::numBCEdges), numCol(mathVar::nGauss1D+1);
        if (Flux)
        {
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfRow "<<numRow<<"\n"<<"NumberOfCol "<<numCol<<"\n"<<"{\n";
            for (int i = 0; i < numRow; i++)
            {
                for (int j = 0; j < numCol; j++)
                {
                    Flux << array[i][j] << " ";
                }
                Flux << "\n";
            }
            Flux<<"}\n";
        }
        else
        {
            std::cout<<"Can not open file "<<loc<<" to write.\n";
            std::cout << "DGSolver will exit after you hit return.\n";
            exit(EXIT_FAILURE);
        }
    }

    void readFromFile(std::string loc, std::string fileName, bool exitWhenFileNotFound, double **array)
    {
        int row(meshVar::numBCEdges), col(1), expectNumCol(mathVar::nGauss1D+1);
        bool isColAvail(true);
        std::ifstream Flux(loc.c_str());
        if (Flux)
        {
            std::string line(" "), checkStr;
            int iElem(0);
            bool startToRead(false);
            //std::tie(row,std::ignore)=auxUlti::lookForDataOfKeyword(loc,"NumberOfRow");
            std::tie(col,isColAvail)=auxUlti::lookForDataOfKeyword(loc,"NumberOfCol");

            if (!isColAvail)
                col=1;

            //Check whether col # nGauss
            bool calcMeanVal(false);
            if (expectNumCol!=col)
                calcMeanVal=true;

            while (std::getline(Flux, line))
            {
                if (startToRead==false)
                {
                    std::istringstream data(line);
                    data >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
                    if (iElem<row)
                    {
                        std::istringstream data(line);
                        if (!calcMeanVal)
                        {
                            for (int i=0; i<col; i++)
                            {
                                data >> array[iElem][i];
                            }
                        }
                        else
                        {
                            double number(0.0),meanVar(0.0),sum(0.0);
                            for (int i=0; i<col; i++)
                            {
                                data >> number;
                                sum+=number;
                            }
                            meanVar=sum/col;

                            for (int j=0; j<expectNumCol; j++)
                            {
                                array[iElem][j] = meanVar;
                            }
                        }
                        iElem++;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
        else
        {
            if (exitWhenFileNotFound)
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(fileName, loc));
            }
        }
    }

    void writeGaussPtsCoor(std::string loc, std::string name)
    {
        std::ofstream Flux(loc.c_str());
        int numRow(meshVar::numBCEdges), numCol(2);
        if (Flux)
        {
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfRow "<<numRow<<"\n"<<"NumberOfCol "<<numCol<<"\n"<<"{\n";
            for (int i = 0; i < numRow; i++)
            {
                Flux << meshVar::GaussPtsOnBCEdge_x[i] << " " << meshVar::GaussPtsOnBCEdge_y[i] << "\n";
            }
            Flux<<"}\n";
        }
        else
        {
            std::cout<<"Can not open file "<<loc<<" to write.\n";
            std::cout << "DGSolver will exit after you hit return.\n";
            exit(EXIT_FAILURE);
        }
    }
}
