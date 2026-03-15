#include "DGIOLib.h"
#include "DGMath.h"
#include "DGMessagesLib.h"
#include "ConstDeclaration.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGProcLib.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib> //include this library to create folder in linux
#include <mpi.h>
#include "DGPostProcessLib.h"
#include "./parallelFunctions/generalParallelFuncs.h"

//Limiters
#include "./limiters/limiterController.h"


//Extended NSF model
#include "./extNSFEqns/FranzDurst/DurstModel.h"
#include "./extNSFEqns/FranzDurst/IO.h"

//Read BC
#include "boundaryConditions/readBCInfor/BCReader.h"

//Non equilibrium BCs
#include "boundaryConditions/customBCs/nonEquilibriumBCs/nonEqmBCs_GenFuncs.h"

namespace IO
{
    /**
     * @brief Function displays DG2DSolver Logo
     */
	void dispLogo()
	{
		std::string logoStr(" ");
		logoStr = message::headerFile();
		std::cout << logoStr << std::endl;
	}

    /**
     * @brief Function gets submitted case's information
     */
	void getCase()
	{
        if (systemVar::currentProc==0)
        {
            std::cout << "***Getting case's information***\n";
        }
		systemVar::wD = auxUlti::workingdir();

		/*Get caseName from SubmitCase*/
		std::string submitLoc(" ");  //submitLoc contents location of submitingCase
        submitLoc = systemVar::wD + "/CASES/SubmitCase.txt";

		std::ifstream submitCaseFlux(submitLoc.c_str());
		if (submitCaseFlux)
		{
			std::string line, keyWord;
            int keyWFlag(0);
			while (std::getline(submitCaseFlux, line))  //Read submitCaseFlux line by line
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					ptr.push_back(keyWord);
				}

                int numWrd = static_cast<int>(ptr.size());

				if (numWrd == 2)
				{
					std::string str1(ptr[0]), str2("caseName");
					if (str1.compare(str2) == 0)
					{
						systemVar::caseName = ptr[1];
						keyWFlag = 1;
                        if (systemVar::currentProc==0)
                        {
                            std::cout << "	Case " << systemVar::caseName << " has been submitted\n";
                        }
					}
				}
			}
			if (keyWFlag == 0)
			{
				std::cout << message::undfKeyW("caseName", submitLoc) << std::endl;
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/"), "", message::opFError("SubmitCase.txt", submitLoc));
		}
        systemVar::pwd = systemVar::wD + "/CASES/" + systemVar::caseName;
	}

    void loadMesh(std::string mode)
	{
        std::string  Elem1DLoc, ptLoc, Elem2DLoc, bcLoc;
		/*Declare loading locations*/
        if (mode.compare("p")!=0)
        {
            Elem1DLoc = systemVar::pwd + "/Constant/Mesh/Elements1D.txt";
            ptLoc = systemVar::pwd + "/Constant/Mesh/Points.txt";
            Elem2DLoc = systemVar::pwd + "/Constant/Mesh/Elements2D.txt";
            bcLoc = systemVar::pwd + "/Constant/boundaryPatch.txt";
        }
        else {
            Elem1DLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/Elements1D.txt";
            ptLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/Points.txt";
            Elem2DLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/Elements2D.txt";
            bcLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/boundaryPatch.txt";
        }

		/*Load Points*/
		std::ifstream ptFlux(ptLoc.c_str());
		if (ptFlux)
		{
            std::string line(" "), checkStr;
			double x, y, z;
            int iPoint(0);
            bool startToRead(false);

            //Doc so luong points (line dau tien cua file)
            std::tie(meshVar::npoin,std::ignore)=auxUlti::lookForDataOfKeyword(ptLoc,"NumberOfEntities");

            //Resize mang meshVar::Points
            meshVar::Points = auxUlti::resize2DArray(meshVar::npoin,3,0.0);
			
            int temp;
			while (std::getline(ptFlux, line))
			{
                //auxUlti::addRowTo2DDoubleArray(meshVar::Points, 3);
                if (startToRead==false)
                {
                    std::istringstream ptData(line);
                    ptData >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
					if (iPoint<meshVar::npoin)
					{
						std::istringstream ptData(line);
                        ptData >> temp >> x >> y >> z;
						meshVar::Points[iPoint][0] = x;
						meshVar::Points[iPoint][1] = y;
						meshVar::Points[iPoint][2] = z;
						iPoint++;
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
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Points.txt", ptLoc));
		}

		/*Load 1D elements*/
		std::ifstream Elem1DFlux(Elem1DLoc.c_str());
		if (Elem1DFlux)
		{
            int node1(0), node2(0), bcGrp(0), i1D(0);  //nelem1D is number of 1D elements
            std::string line(" "), checkStr(" ");
            bool startToRead(false);

            //Doc so luong element 1D
            std::tie(meshVar::nelem1D,std::ignore)=auxUlti::lookForDataOfKeyword(Elem1DLoc,"NumberOfEntities");

            //Resize mang meshVar::Elements1D
            meshVar::Elements1D = auxUlti::resize2DIntArray(meshVar::nelem1D,3,-1);

            int temp;
			while (std::getline(Elem1DFlux, line))
			{
                //auxUlti::addRowTo2DIntArray(meshVar::Elements1D, 3);
                if (startToRead==false)
                {
                    std::istringstream elData(line);
                    elData >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
					if (i1D<meshVar::nelem1D)
					{
						std::istringstream elData(line);
                        elData >> temp>> node1 >> node2 >> bcGrp;
						meshVar::Elements1D[i1D][0] = node1 - 1;
						meshVar::Elements1D[i1D][1] = node2 - 1;
						meshVar::Elements1D[i1D][2] = bcGrp;
						i1D++;
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
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Elements1D.txt", Elem1DLoc));
		}

		/*Load 2D elements*/
		std::ifstream Elem2DFlux(Elem2DLoc.c_str());
		if (Elem2DFlux)
		{
            int node1(0), node2(0), node3(0), node4(0), i2D(0);  //nelem2D is number of 2D elements
            std::string line(" "), checkStr;
            bool startToRead(false);

            //Doc so luong element 2D
            std::tie(meshVar::nelem2D,std::ignore)=auxUlti::lookForDataOfKeyword(Elem2DLoc,"NumberOfEntities");

            //Resize mang meshVar::Elements2D
            meshVar::Elements2D = auxUlti::resize2DIntArray(meshVar::nelem2D,4,-1);

            int temp;
			while (std::getline(Elem2DFlux, line))
			{
                //auxUlti::addRowTo2DIntArray(meshVar::Elements2D, 4);
                if (startToRead==false)
                {
                    std::istringstream elData(line);
                    elData >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
					if (i2D<meshVar::nelem2D)
					{
						std::istringstream elData(line);
                        elData >> temp>> node1 >> node2 >> node3 >> node4;
						meshVar::Elements2D[i2D][0] = node1 - 1;
						meshVar::Elements2D[i2D][1] = node2 - 1;
						meshVar::Elements2D[i2D][2] = node3 - 1;
						meshVar::Elements2D[i2D][3] = node4 - 1;
						i2D++;
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
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Elements2D.txt", Elem2DLoc));
		}

        /*Load boundary conditions from boundaryPatch file*/
		std::ifstream bcFlux(bcLoc.c_str());
		if (bcFlux)
		{
            int boundIndex(0);
			std::string line(" "), keyWord(" ");
            //Doc so luong boundary
            std::tie(meshVar::nBc,std::ignore)=auxUlti::lookForDataOfKeyword(bcLoc,"NumberOfEntities");
			
			//bien nBc se tang them 1 (vi them 1 bien matched) khi readingMode parallel
            if (mode.compare("p")==0)
            {
                meshVar::nBc++;
            }
			
            meshVar::BoundaryType = auxUlti::resize2DIntArray(meshVar::nBc,3,0);

			while (std::getline(bcFlux, line))
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					if ((keyWord.compare(" ") != 0))
					{
						ptr.push_back(keyWord);
					}
				}

                int numWrd = static_cast<int>(ptr.size());
				if (numWrd != 0)
				{
					std::string str1(ptr[0]);
					if (str1.compare("{") == 0)  //0 means str1 is the same as str2
					{
						//start group
						boundIndex++;
					}

					if ((numWrd >= 2)&(boundIndex != 0))
					{
						if ((boundIndex <= bcSize))
						{
							meshVar::BoundaryType[boundIndex - 1][0] = boundIndex;
							std::string str2 = ptr[1];
							if (str1.compare("Type") == 0)
							{
								if (str2.compare("wall") == 0)
								{
                                    meshVar::BoundaryType[boundIndex - 1][1] = meshVar::BCTypeID::wall;  //type 1 is wall
								}
								else if (str2.compare("patch") == 0)
								{
                                    meshVar::BoundaryType[boundIndex - 1][1] = meshVar::BCTypeID::patch;  //type 2 is patch
								}
								else if (str2.compare("symmetry") == 0)
								{
                                    meshVar::BoundaryType[boundIndex - 1][1] = meshVar::BCTypeID::symmetry;  //type 3 is symmetry
								}
                                else if (str2.compare("matched") == 0)
                                {
                                    meshVar::BoundaryType[boundIndex - 1][1] = meshVar::BCTypeID::matched;  //type 4 is matched (boundary condition type for parallel computing)
                                }
								else
								{
									std::string errorStr = "Boundary type <" + ptr[1] + R"(> is unknown, available boundary types in file boundaryPatch are:
	wall
	patch
    symmetry
    matched)";
									message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
								}
							}
						}
						else
						{
                            std::string errorStr = "Maximum number of boundaries DG2D solver supports is only " + std::to_string(bcSize) + ". Make sure the number of boundaries is less than " + std::to_string(bcSize);
							message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
						}
					}
				}
                //meshVar::nBc = boundIndex;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("boundaryPatch.txt", bcLoc));
		}

        if (mode.compare("p")==0)
        {
            int meshConnectionLength(0);
            bool startToRead(false);
            std::string mshConnect = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/meshConnection.txt",
                    checkStr;
            std::ifstream mshConnectFlux(mshConnect.c_str());
            if (mshConnectFlux)
            {
                std::string line("");
                int counter = 0;
                std::tie(meshConnectionLength,std::ignore)=auxUlti::lookForDataOfKeyword(mshConnect,"NumberOfEntities");
                meshVar::meshConnection = auxUlti::resize2DIntArray(meshConnectionLength,3,0);

                while (std::getline(mshConnectFlux, line))
                {
                    //auxUlti::addRowTo2DIntArray(meshVar::meshConnection, 3);
                    if (startToRead==false)
                    {
                        std::istringstream Data(line);
                        Data >> checkStr;
                        if (checkStr.compare("{") == 0)
                        {
                            startToRead=true;
                        }
                    }
                    else
                    {
						if (counter<meshConnectionLength)
						{
							std::istringstream Data(line);
							Data >> meshVar::meshConnection[counter][0] >> meshVar::meshConnection[counter][1] >> meshVar::meshConnection[counter][2];
							counter++;
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
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("meshConnection.txt", mshConnect));
            }
        }

        //Read sendRecvOrder
        if (mode.compare("p")==0)
        {
            bool startToRead(false);
            std::string sendRecvOrderLoc = systemVar::pwd + "/Constant/Mesh/sendRecvOrder.txt",
                    checkStr;
            std::ifstream sendRecvOrderFlux(sendRecvOrderLoc.c_str());
            if (sendRecvOrderFlux)
            {
                std::string line("");
                int counter = 0;
                std::tie(systemVar::sendRecvOrder_length,std::ignore)=auxUlti::lookForDataOfKeyword(sendRecvOrderLoc,"NumberOfEntities");
                systemVar::sendRecvOrder = auxUlti::resize2DIntArray(systemVar::sendRecvOrder_length,2,0);

                while (std::getline(sendRecvOrderFlux, line))
                {
                    //auxUlti::addRowTo2DIntArray(meshVar::meshConnection, 3);
                    if (startToRead==false)
                    {
                        std::istringstream Data(line);
                        Data >> checkStr;
                        if (checkStr.compare("{") == 0)
                        {
                            startToRead=true;
                        }
                    }
                    else
                    {
                        if (counter<systemVar::sendRecvOrder_length)
                        {
                            std::istringstream Data(line);
                            Data >> systemVar::sendRecvOrder[counter][0] >> systemVar::sendRecvOrder[counter][1];
                            counter++;
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
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("sendRecvOrder.txt", sendRecvOrderLoc));
            }
        }
	}

    void SaveMeshInfor(std::string mode)
	{
		/*List of arrays need to be saved:
		- inedel: array contents information of edges belong to elements, column index is element index, each row in column contents index of edge belong to element, number of row is 4 because of default quad element.
		- ineled: array contents information of elements sharing edges, column index is edge index, each row in column contents index of element which edge is belong to, row 3 contents pointer.
		- inpoed: array contents information of points belong to edges, column 3 contents group which edge belongs to (group 0 is internal group), column 4 contents type of boundary (type 0 is internal edge)
		- normalVector: array contents information of normal vector of edges
		- MasterElemOfEdge: array content master element of iedge, use it with normalVector to get information of normal vector of edge*/

        std::string inedelLoc, ineledLoc, inpoedLoc, normVectorLoc, MasterElemOfEdgeLoc, Elem1DLoc, ptLoc, Elem2DLoc, bcLoc;
		
        /*Declare loading locations*/
        if (mode.compare("p")!=0)
        {
            inedelLoc = systemVar::pwd + "/Constant/Mesh/inedel.txt";
            ineledLoc = systemVar::pwd + "/Constant/Mesh/ineled.txt";
            inpoedLoc = systemVar::pwd + "/Constant/Mesh/inpoed.txt";
            normVectorLoc = systemVar::pwd + "/Constant/Mesh/normalVector.txt";
            MasterElemOfEdgeLoc = systemVar::pwd + "/Constant/Mesh/MasterElemOfEdge.txt";
        }
        else {
            inedelLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/inedel.txt";
            ineledLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/ineled.txt";
            inpoedLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/inpoed.txt";
            normVectorLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/normalVector.txt";
            MasterElemOfEdgeLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/MasterElemOfEdge.txt";
        }

		/*inedel*/
		std::ofstream Fluxinedel(inedelLoc.c_str());
		if (Fluxinedel)
		{
			Fluxinedel << systemVar::headerFile << std::endl << "	inedel array\n" << "\n";
            for (int i = 0; i < meshVar::nelem2D; i++)
			{
                for (int j = 0; j < 4; j++)
				{
					Fluxinedel << meshVar::inedel[i][j] << " ";
				}
				Fluxinedel << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("inedel.txt", inedelLoc));
		}

		/*ineled*/
		std::ofstream Fluxineled(ineledLoc.c_str());
		if (Fluxineled)
		{
			Fluxineled << systemVar::headerFile << std::endl << "	ineled array\n" << "\n";
            for (int i = 0; i < meshVar::inpoedCount; i++)
			{
                for (int j = 0; j < 2; j++)
				{
					Fluxineled << meshVar::ineled[i][j] << " ";
				}
				Fluxineled << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("ineled.txt", ineledLoc));
		}

		/*inpoed*/
		std::ofstream Fluxinpoed(inpoedLoc.c_str());
		if (Fluxinpoed)
		{
			Fluxinpoed << systemVar::headerFile << std::endl << "	inpoed array\n" << "\n";
            for (int i = 0; i < meshVar::inpoedCount; i++)
			{
                for (int j = 0; j < 5; j++)
				{
					Fluxinpoed << meshVar::inpoed[i][j] << " ";
				}
				Fluxinpoed << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("inpoed.txt", inpoedLoc));
		}

		/*normalVector*/
		std::ofstream FluxnormVector(normVectorLoc.c_str());
		if (FluxnormVector)
		{
			FluxnormVector << systemVar::headerFile << std::endl << "	normalVector array\n" << "\n";
			for (int i = 0; i < meshVar::inpoedCount; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					FluxnormVector << meshVar::normalVector[i][j] << " ";
				}
				FluxnormVector << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("normalVector.txt", normVectorLoc));
		}

		/*inpoed*/
		std::ofstream FluxMasterElemOfEdge(MasterElemOfEdgeLoc.c_str());
		if (FluxMasterElemOfEdge)
		{
			FluxMasterElemOfEdge << systemVar::headerFile << std::endl << "	MasterElemOfEdge array\n" << "\n";
			for (int i = 0; i < meshVar::inpoedCount; i++)
			{
				FluxMasterElemOfEdge << meshVar::MasterElemOfEdge[i] << " ";
				FluxMasterElemOfEdge << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("MasterElemOfEdge.txt", MasterElemOfEdgeLoc));
		}
	}

    /**
     * @brief Function reads number of assigned cores for running parallel
     */
    void readNumberOfCores()
    {
        /*Read totalProcess in DGOptions*/
        std::string DGOptfileName("DGOptions.txt");
        std::string DGOptLoc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
        std::string DGOptkeyWordsDouble[1] = {" "}, DGOptkeyWordsInt[1] = {"totalProcess"}, DGOptkeyWordsBool[1] = {" "}, DGOptkeyWordsStr[1] = {" "};
        double DGOptoutDB[1] = {};
        int DGOptoutInt[1] = {};
        bool DGOptoutBool[1] = {};
        std::string DGOptoutStr[1] = {};

        readDataFile(DGOptfileName, DGOptLoc, DGOptkeyWordsDouble, DGOptkeyWordsInt, DGOptkeyWordsBool, DGOptkeyWordsStr, DGOptoutDB, DGOptoutInt, DGOptoutBool, DGOptoutStr, 0, 1, 0, 0);

        systemVar::totalProc = DGOptoutInt[0];
    }



	void readDataFile(std::string fileName, std::string direction, std::string keyWordsDbl[], std::string keyWordsInt[], std::string keyWordsBool[], std::string keyWordsStr[], double *outDbl, int *outInt, bool *outBool, std::string *outStr, int numParamDbl, int numParamInt, int numParamBool, int numParamStr)  //Declaration of funciton which returns pointer of 1D array
	{
		/*User's guide:
		This function returns datas of type double and type int read from files.
		Input arguments:
		- fileName: name of file you want to read (file name and extension)
        - direction: working diectory contents data file (with no "/" characters at the end)
		- keyWordsDbl: array contents keyWords of double values listed in file
		- keyWordsInt: array contents keyWords of int values listed in file
		- keyWordsBool: array contents keyWords of bool values listed in file
		- keyWordsStr: array contents keyWords of string values listed in file
		- outDbl: output array contents double values
		- outInt: output array contents int values
        - numParamDbl, numParamInt, numParamBool, numParamStr: number of double, int, bool, string parameters

        To write comment, use the following syntax
        // comment... */

        //Modify number of parameters
        if (numParamDbl<0) numParamDbl=0;
        if (numParamInt<0) numParamInt=0;
        if (numParamBool<0) numParamBool=0;
        if (numParamStr<0) numParamStr=0;

        //Checking array
        std::vector<bool> check_double(numParamDbl,false),
                check_int(numParamInt,false),
                check_bool(numParamBool,false),
                check_str(numParamStr,false);

		double dataDbl(0.0);
        int dataInt(0);
		std::string dataStr("abc");
        /*
        if (systemVar::currentProc==0)
        {
            std::cout << "	Reading " << fileName <<"\n";
        }*/

        std::string FileLoc(direction + "/" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
		int indexDbl(0), indexInt(0), indexBool(0), indexStr(0);
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
					std::string str1(ptr[0]), str2(" "), str3(" "), str4(" "), str5(" ");
					std::istringstream strdata(ptr[1]);

                    //Doc rieng tung loai data
                    //double
                    if (indexDbl<numParamDbl)
                    {
                        for (int i=0; i<numParamDbl; i++)
                        {
                            if (!check_double[i])
                            {
                                str2 = keyWordsDbl[i];
                                if (str1.compare(str2) == 0)
                                {
                                    strdata >> dataDbl;
                                    outDbl[i]=dataDbl;
                                    indexDbl++;
                                    check_double[i]=true;
                                    goto label;
                                }
                            }
                        }
                    }

                    //int
                    if (indexInt<numParamInt)
                    {
                        for (int i=0; i<numParamInt; i++)
                        {
                            if (!check_int[i])
                            {
                                str2 = keyWordsInt[i];
                                if (str1.compare(str2) == 0)
                                {
                                    strdata >> dataInt;
                                    outInt[i]=dataInt;
                                    indexInt++;
                                    check_int[i]=true;
                                    goto label;
                                }
                            }
                        }
                    }

                    //string
                    if (indexStr<numParamStr)
                    {
                        for (int i=0; i<numParamStr; i++)
                        {
                            if (!check_str[i])
                            {
                                str2 = keyWordsStr[i];
                                if (str1.compare(str2) == 0)
                                {
                                    strdata >> dataStr;
                                    outStr[i]=dataStr;
                                    indexStr++;
                                    check_str[i]=true;
                                    goto label;
                                }
                            }
                        }
                    }

                    //bool
                    if (indexBool<numParamBool)
                    {
                        for (int i=0; i<numParamBool; i++)
                        {
                            if (!check_bool[i])
                            {
                                str2 = keyWordsBool[i];
                                if (str1.compare(str2) == 0)
                                {
                                    if ((ptr[1].compare("true") == 0) || (ptr[1].compare("yes") == 0))
                                    {
                                        outBool[i] = true;
                                    }
                                    else if ((ptr[1].compare("false") == 0) || (ptr[1].compare("no") == 0))
                                    {
                                        outBool[i] = false;
                                    }
                                    indexBool++;
                                    check_bool[i]=true;
                                    goto label;
                                }
                            }
                        }
                    }
                }
                label:
				if ((indexDbl>=numParamDbl) && (indexInt>=numParamInt) && (indexBool >= numParamBool) && (indexStr >= numParamStr))
				{
					break;
				}
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}

        //Check whether reading request is fulfill
        for (int i=0; i<numParamDbl; i++)
        {
            if (!check_double[i])
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::undfKeyW(keyWordsDbl[i], fileName));
            }
        }
        for (int i=0; i<numParamInt; i++)
        {
            if (!check_int[i])
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::undfKeyW(keyWordsInt[i], fileName));
            }
        }
        for (int i=0; i<numParamBool; i++)
        {
            if (!check_bool[i])
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::undfKeyW(keyWordsBool[i], fileName));
            }
        }
        for (int i=0; i<numParamStr; i++)
        {
            if (!check_str[i])
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::undfKeyW(keyWordsStr[i], fileName));
            }
        }
	}

    std::tuple<bool, double**, int> read2DArray(int column, std::string location, std::string fileName, bool exitWhenFileNotFound)
    {
        int length(0);
        double**array;
        bool isFileAvailable;
        std::ifstream Flux(location.c_str());
        if (Flux)
        {
            isFileAvailable=true;
            std::string line(" "), checkStr;
            int iElem(0);
            bool startToRead(false);
            std::tie(length,std::ignore)=auxUlti::lookForDataOfKeyword(location,"NumberOfEntities");
            array=auxUlti::resize2DArray(length,column,0.0);

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
                    if (iElem<length)
                    {
                        std::istringstream data(line);
                        for (int i=0; i<column; i++)
                        {
                            data >> array[iElem][i];
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
            isFileAvailable=false;
            if (exitWhenFileNotFound)
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(fileName, location));
            }
        }
        return std::make_tuple(isFileAvailable,array,length);
    }

    std::tuple<bool, int**, int> read2DIntArray(int column, std::string location, std::string fileName, bool exitWhenFileNotFound)
    {
        int length(0);
        int**array;
        bool isFileAvailable;
        std::ifstream Flux(location.c_str());
        if (Flux)
        {
            isFileAvailable=true;

            std::string line(" "), checkStr;
            int iElem(0);
            bool startToRead(false);
            std::tie(length,std::ignore)=auxUlti::lookForDataOfKeyword(location,"NumberOfEntities");
            array=auxUlti::resize2DIntArray(length,column,0);

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
                    if (iElem<length)
                    {
                        std::istringstream data(line);
                        for (int i=0; i<column; i++)
                        {
                            data >> array[iElem][i];
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
            isFileAvailable=false;
            if (exitWhenFileNotFound)
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(fileName, location));
            }
        }
        return std::make_tuple(isFileAvailable,array,length);
    }

    std::tuple<bool, double*, int> read1DArray(std::string location, std::string fileName, bool exitWhenFileNotFound)
    {
        int length(0);
        double*array;
        bool isFileAvailable;
        std::ifstream Flux(location.c_str());
        if (Flux)
        {
            isFileAvailable=true;

            std::string line(" "), checkStr;
            int iElem(0);
            bool startToRead(false);
            std::tie(length,std::ignore)=auxUlti::lookForDataOfKeyword(location,"NumberOfEntities");
            array = new double [length];

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
                    if (iElem<length)
                    {
                        std::istringstream data(line);
                        data >> array[iElem];
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
            isFileAvailable=false;
            if (exitWhenFileNotFound)
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(fileName, location));
            }
        }
        return std::make_tuple(isFileAvailable,array,length);
    }

    std::tuple<bool,int*, int> read1DIntArray(std::string location, std::string fileName, bool exitWhenFileNotFound)
    {
        int length(0);
        int*array;
        bool isFileAvailable;
        std::ifstream Flux(location.c_str());
        if (Flux)
        {
            isFileAvailable=true;

            std::string line(" "), checkStr;
            int iElem(0);
            bool startToRead(false);
            std::tie(length,std::ignore)=auxUlti::lookForDataOfKeyword(location,"NumberOfEntities");
            array = new int [length];

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
                    if (iElem<length)
                    {
                        std::istringstream data(line);
                        data >> array[iElem];
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
            isFileAvailable=false;
            if (exitWhenFileNotFound)
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(fileName, location));
            }
        }
        return std::make_tuple(isFileAvailable,array,length);
    }

	void residualOutput(double rhoResGlobal, double rhouResGlobal, double rhovResGlobal, double rhoEResGlobal)
	{
        if (systemVar::iterCount<5)
		{
			systemVar::rhoResNorm = rhoResGlobal;
			systemVar::rhouResNorm = rhouResGlobal;
			systemVar::rhovResNorm = rhovResGlobal;
			systemVar::rhoEResNorm = rhoEResGlobal;
			rhoResGlobal = 1.0;
			rhouResGlobal = 1.0;
			rhovResGlobal = 1.0;
			rhoEResGlobal = 1.0;
		}
		else
		{
			rhoResGlobal /= systemVar::rhoResNorm;
			rhouResGlobal /= systemVar::rhouResNorm;
			rhovResGlobal /= systemVar::rhovResNorm;
			rhoEResGlobal /= systemVar::rhoEResNorm;
		}

        std::cout << "Time step: " << dt << std::endl;
        std::cout << "Residuals: ddt(rho)=" << rhoResGlobal << ", ddt(rhou)=" << rhouResGlobal << ", ddt(rhov)=" << rhovResGlobal << ", ddt(rhoE)=" << rhoEResGlobal << std::endl;
        std::cout << std::endl;

        IO::writeResiduals(systemVar::iterCount,rhoResGlobal,rhouResGlobal,rhovResGlobal,rhoEResGlobal);
	}

    void readDecomposedMeshInfor()
    {
        std::string rankOf2DElemLoc(systemVar::pwd+"/Constant/Mesh/rankOf2DElem.txt"),
                Elem2DlocalIdWithRankLoc(systemVar::pwd+"/Constant/Mesh/Elem2DlocalIdWithRank.txt");

        std::ifstream File1Flux(rankOf2DElemLoc.c_str());
        int i(0);
        if (File1Flux)
        {
            File1Flux>>meshVar::rankOf2DElem[i];
            i++;
        }

        i=0;
        std::ifstream File2Flux(Elem2DlocalIdWithRankLoc.c_str());
        if (File2Flux)
        {
            File2Flux>>meshVar::Elem2DlocalIdWithRank[i];
            i++;
        }
    }

    //Ham write cac field roi rac (discreted) xuong file fileName tai folder Loc
    void writeDiscretedFields(std::string Loc, std::string fileName, double **Var)
    {
        std::string fileLoc(Loc + "/" + fileName);
        std::ofstream fileFlux(fileLoc.c_str());
        fileFlux<<systemVar::headerFile<<"NumberOfEntities "<<meshVar::nelem2D<<"\n"<<"{\n";
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
        {
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
            {
                fileFlux << Var[nelem][iorder] << " ";
            }
            fileFlux << "\n";
        }
        fileFlux<<"}";
    }

    //Ham read cac field roi rac (discreted) xuong file fileName tai folder Loc
    void readDiscretedFields(std::string Loc, std::string fileName, double **Var)
    {
        std::string fileLoc(Loc + "/" + fileName), checkStr;
        std::ifstream FileFlux(fileLoc.c_str());

        bool startToRead(false);
        if (FileFlux)
        {
            int nelement(0);
            std::string line;
            while (std::getline(FileFlux, line))
            {
                //auxUlti::addRowTo2DIntArray(meshVar::Elements1D, 3);
                if (startToRead==false)
                {
                    std::istringstream lineflux(line);
                    lineflux >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
                    if (nelement<meshVar::nelem2D)
                    {
                        std::istringstream lineflux(line);
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            lineflux >> Var[nelement][iorder];
                        }
                        nelement++;
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
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
        }
    }

    void saveCase()
	{
		std::string iter_str = std::to_string(systemVar::iterCount);
        std::string fileName, Loc, fileLoc;
        Loc = auxUlti::createTimeStepFolder(systemVar::iterCount,"case");

		/*Conservative variables*/
        IO::writeDiscretedFields(Loc,"rho.txt",rho);
        IO::writeDiscretedFields(Loc,"rhou.txt",rhou);
        IO::writeDiscretedFields(Loc,"rhov.txt",rhov);
        IO::writeDiscretedFields(Loc,"rhoE.txt",rhoE);
		/*end of saving conservative variables*/

        //Non Equilibrium: write surface values
        nonEquilibriumBCs_IO::writeSurfaceValues(Loc);

        //Write variables on wall
        postProcessing_Surface::writeVarsAtWall(Loc);

		/*Residual normalized coeffs*/
		fileName = "ResidualNormCoeffs.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ofstream fileFluxResNorm(fileLoc.c_str());
        fileFluxResNorm << systemVar::rhoResNorm << " " << systemVar::rhouResNorm << " " << systemVar::rhovResNorm << " " << systemVar::rhoEResNorm << "\n";

		/*Time informations*/
        if (systemVar::currentProc==0)
        {
            Loc = systemVar::wD + "/CASES/" + systemVar::caseName;
            fileName = "time.txt";
            fileLoc = (Loc + "/" + fileName);
            std::ofstream fileFluxTime(fileLoc.c_str());
            fileFluxTime << systemVar::iterCount << " "<<runTime<< "\n";
        }

        MPI_Barrier(MPI_COMM_WORLD);
	}

    void saveCase_reconstruct()
    {
        std::string iter_str = std::to_string(systemVar::iterCount);
        std::string fileName, Loc, fileLoc;
        Loc = auxUlti::createTimeStepFolder(systemVar::iterCount,"case");

        /*Conservative variables*/
        IO::writeDiscretedFields(Loc,"rho.txt",rho);
        IO::writeDiscretedFields(Loc,"rhou.txt",rhou);
        IO::writeDiscretedFields(Loc,"rhov.txt",rhov);
        IO::writeDiscretedFields(Loc,"rhoE.txt",rhoE);
        /*end of saving conservative variables*/
    }

    void loadCase(std::string mode)
	{
        //Read file time.txt
        std::string fileName("time.txt"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName);
        std::string fileLoc(Loc + "/" + fileName);
		std::ifstream FileFluxTime(fileLoc.c_str());
		if (FileFluxTime)
		{
			std::string line;
			std::getline(FileFluxTime, line);
			std::istringstream lineflux(line);
            lineflux >> systemVar::iterCount>>runTime;
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}

		std::string iter_str = std::to_string(systemVar::iterCount);
        Loc = systemVar::wD + "/CASES/" + systemVar::caseName;
        if (mode.compare("p")==0)
        {
            Loc=Loc+"/Processor"+std::to_string(systemVar::currentProc) + "/" + iter_str;
        }
        else {
            Loc=Loc+"/" +iter_str;
        }

        IO::readDiscretedFields(Loc,"rho.txt",rho);
        IO::readDiscretedFields(Loc,"rhou.txt",rhou);
        IO::readDiscretedFields(Loc,"rhov.txt",rhov);
        IO::readDiscretedFields(Loc,"rhoE.txt",rhoE);

        //Non Equilibrium: read surface values
        nonEquilibriumBCs_IO::readSurfaceValues(Loc);

		//Read residual norm coeffs
		fileName = "ResidualNormCoeffs.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ifstream FileFluxResNorm(fileLoc.c_str());
		if (FileFluxResNorm)
		{
			std::string line;
			std::getline(FileFluxResNorm, line);
			std::istringstream lineflux(line);
			lineflux >> systemVar::rhoResNorm >> systemVar::rhouResNorm >> systemVar::rhovResNorm >> systemVar::rhoEResNorm;
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}
	}

    void loadpTU(std::string mode)
    {
        /*Read U*/  //U must be read first of all
        readVectorBC::u(mode);

        /*Read p*/
        readScalarBC::p(mode);

        /*Read T*/
        readScalarBC::T(mode);

        /*Tinh so Mach*/
        flowProperties::Mach=pow(pow(iniValues::uIni,2)+pow(iniValues::vIni,2)+pow(iniValues::wIni,2),0.5)/math::CalcSpeedOfSound(iniValues::TIni);
    }

    void loadTime()
    {
        //Read file time.txt
        std::string fileName("time.txt"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName);
        std::string fileLoc(Loc + "/" + fileName);
        std::ifstream FileFluxTime(fileLoc.c_str());
        if (FileFluxTime)
        {
            std::string line;
            std::getline(FileFluxTime, line);
            std::istringstream lineflux(line);
            lineflux >> systemVar::iterCount>>runTime;
        }
        else
        {
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
        }
    }

    void write2DDoubleArrayToFile_typeVector(std::vector<std::vector<double>> &array, std::string loc, std::string name, int numRow, int numCol)
    {
        std::ofstream Flux(loc.c_str());
        if (Flux)
        {
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfEntities "<<numRow<<"\n"<<"{\n";
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

    void write2DDoubleArrayToFile_typeArray(double **array, std::string loc, std::string name, int numRow, int numCol)
    {
        std::ofstream Flux(loc.c_str());
        if (Flux)
        {
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfEntities "<<numRow<<"\n"<<"{\n";
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

    void writeResiduals(int iter, double rhoRes, double rhouRes, double rhovRes, double rhoERes)
    {
        IO::openFileToAppend((systemVar::pwd+"/ResidualOutput.txt"),std::to_string(iter)+" "
                             +std::to_string(rhoRes)+" "
                             +std::to_string(rhouRes)+" "
                             +std::to_string(rhovRes)+" "
                             +std::to_string(rhoERes)+"\n");
    }

    void write2DIntArrayToFile_typeVector(std::vector<std::vector<int>> &array, std::string loc, std::string name, int numRow, int numCol)
	{
		std::ofstream Flux(loc.c_str());
		if (Flux)
		{
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfEntities "<<numRow<<"\n"<<"{\n";
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

    void write2DIntArrayToFile_typeArray(int **array, std::string loc, std::string name, int numRow, int numCol)
    {
        std::ofstream Flux(loc.c_str());
        if (Flux)
        {
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfEntities "<<numRow<<"\n"<<"{\n";
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

    void openFileToAppend(std::string Loc, std::string content)
    {
        std::fstream fs;
          fs.open (Loc, std::fstream::app);

          fs << content;

          fs.close();
    }

	namespace importCase
	{
    //Phan nay k can
		void importResultsFromAnotherCase()
		{
            //SUA TAM THOI
            //import case parallel
            std::string sourceCaseName(" "), sourceLoc(" "), timeLoc(" "), fileName(" "), fileLoc(" "), processor("");

            if (systemVar::currentProc==0)
            {
                std::cout << "NOTE: source case must has the same mesh with current case.\n" << "Enter source name: ";
                std::cin>>sourceCaseName;
                //Send command to all processors
                for (int irank=1;irank<systemVar::totalProc;irank++) {
                    parallelFuncs_Gen::sendString(sourceCaseName,irank,2);
                }
            }
            else {
                sourceCaseName=parallelFuncs_Gen::receiveString(0,2);
            }

            sourceLoc = systemVar::wD + "/CASES/" + sourceCaseName;
            timeLoc = sourceLoc + "/time.txt";

            //read source DGOptions
            /*Read DGOptions*/
            std::string DGOptfileName("DGOptions.txt");
            std::string DGOptLoc(systemVar::wD + "/CASES/" + sourceCaseName + "/System");
            std::string DGOptkeyWordsDouble[1] = {}, DGOptkeyWordsInt[1] = {}, DGOptkeyWordsBool[1] = {}, DGOptkeyWordsStr[1] = {"runningMode"};
            double DGOptoutDB[1] = {};
            int DGOptoutInt[1] = {};
            bool DGOptoutBool[1] = {};
            std::string DGOptoutStr[1] = {};

            readDataFile(DGOptfileName, DGOptLoc, DGOptkeyWordsDouble, DGOptkeyWordsInt, DGOptkeyWordsBool, DGOptkeyWordsStr, DGOptoutDB, DGOptoutInt, DGOptoutBool, DGOptoutStr, 0, 0, 0, 1);
            if ((DGOptoutStr[0].compare("parallel") == 0)&& systemVar::parallelMode)
            {
                processor="Processor"+std::to_string(systemVar::currentProc);
            }
            else if (DGOptoutStr[0].compare("sequence") == 0 && !systemVar::parallelMode)
            {
                processor="";
            }
            else {
                std::cout << "runningMode option of source and destination cases are not the same.\n";
                std::cout << "DGSolver will exit after you hit return.\n";
                exit(EXIT_FAILURE);
            }

			//get time at source folder
			std::ifstream timeFlux(timeLoc.c_str());
			if (timeFlux)
			{
				int time(0);
				std::string line;
				std::getline(timeFlux, line);
				std::istringstream lineflux(line);
				lineflux >> time;

                if (systemVar::currentProc==0)
                {
                    std::cout << "Source time: " << std::to_string(time) <<std::endl;
                    std::cout << "Mapping fields:\n";
                }
				//read results from source folder
				fileName = "rho.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
                //mappSourceToCurrent(fileLoc, rho);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rho\n";
                }

				fileName = "rhou.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
                //mappSourceToCurrent(fileLoc, rhou);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhou\n";
                }

				fileName = "rhov.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
                //mappSourceToCurrent(fileLoc, rhov);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhov\n";
                }

				fileName = "rhoE.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
                //mappSourceToCurrent(fileLoc, rhoE);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhoE\n";
                    std::cout << "DONE\n";
                }
			}
			else
			{
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError("time.txt", sourceLoc));
			}
		}

		void mappSourceToCurrent(std::string fileLoc, std::vector<std::vector<double>> &currentResult)
		{
			std::ifstream FileFluxRho(fileLoc.c_str());
			std::vector<double> iniValue(mathVar::orderElem + 1, 0.0);
			if (FileFluxRho)
			{
				int nelement(0);
				double value(0.0);
				std::string line;
				while (std::getline(FileFluxRho, line))
				{
					std::istringstream lineflux(line);
					lineflux >> value;
					iniValue = process::calcIniValues(value, nelement);
					for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
					{
						currentResult[nelement][iorder] = iniValue[iorder];
					}
					nelement++;
				}
			}
			else
			{
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileLoc, fileLoc));
			}
		}
	}

    void write2DDoubleVectorToFile(std::string location, std::string fileName, std::vector<std::vector<double>> &vector)
    {
        std::ofstream FileFlux((location+"/"+fileName).c_str());
        if (FileFlux)
        {
            FileFlux<<message::headerFile()<<fileName<<"\n"<<"NumberOfEntities "<<vector.size()<<"\n"<<"{\n";
            for (int irow=0; irow<static_cast<int>(vector.size()); irow++)
            {
                for (int icol=0; icol<static_cast<int>(vector[irow].size()); icol++)
                {
                    FileFlux<<vector[irow][icol]<<" ";
                }
                FileFlux<<"\n";
            }
            FileFlux<<"}\n";
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(location+"/"+fileName)<<" to write.\n";
        }
    }

    void write1DDoubleVectorToFile(std::string location, std::string fileName, double *vector, int length)
    {
        std::ofstream FileFlux((location+"/"+fileName).c_str());
        FileFlux<<message::headerFile()<<fileName<<"\n"<<"NumberOfEntities "<<length<<"\n"<<"{\n";
        if (FileFlux)
        {
            for (int irow=0; irow<length; irow++)
            {
                FileFlux<<vector[irow]<<"\n";
            }
            FileFlux<<"}\n";
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(location+"/"+fileName)<<" to write.\n";
        }
    }

    void write1DIntVectorToFile(std::string location, std::string fileName, int *vector, int length)
    {
        std::ofstream FileFlux((location+"/"+fileName).c_str());
        if (FileFlux)
        {
            FileFlux<<message::headerFile()<<fileName<<"\n"<<"NumberOfEntities "<<length<<"\n"<<"{\n";
            for (int irow=0; irow<length; irow++)
            {
                FileFlux<<vector[irow]<<"\n";
            }
            FileFlux<<"}\n";
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(location+"/"+fileName)<<" to write.\n";
        }
    }

    void write1DBoolVectorToFile(std::string location, std::string fileName, bool *vector, int length)
    {
        std::ofstream FileFlux((location+"/"+fileName).c_str());
        if (FileFlux)
        {
            FileFlux<<message::headerFile()<<fileName<<"\n"<<"NumberOfEntities "<<length<<"\n"<<"{\n";
            for (int irow=0; irow<length; irow++)
            {
                FileFlux<<vector[irow]<<"\n";
            }
            FileFlux<<"}\n";
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(location+"/"+fileName)<<" to write.\n";
        }
    }

    namespace loadSettingFiles
    {
        /**
         * @brief Function reads constants at initial time
         *
         * - DGOptions: case settings
         * - Material: fluid properties
         * - DGSchemes: numerical schemes
         * - FlowProperties: flow condition: inviscid, viscous, viscous and self-diffusion
         * - Limiters: limiter setting
         * - TBounds: Temperature limit
         */
        void loadConstants()
        {
            IO::loadSettingFiles::DGOptions();
            IO::loadSettingFiles::Material();
            IO::loadSettingFiles::DGSchemes();
            IO::loadSettingFiles::FlowProperties();
            limiter::IOLimiter::readSelectedLimiters();
            IO::loadSettingFiles::TBounds();
        }

        /**
         * @brief Function reads constants while running
         *
         * - DGOptions: case settings
         * - Material: fluid properties
         * - DGSchemes: numerical schemes
         */
        void loadConstantsWhileRunning()
        {
            IO::loadSettingFiles::DGOptions();
            IO::loadSettingFiles::DGSchemes();
            IO::loadSettingFiles::Material();
            //limiter::IOLimiter::readSelectedLimiters();
            //IO::loadSettingFiles::TBounds();
        }

        /**
         * @brief Function reads DGOptions file
         */
        void DGOptions()
        {
            /*Read DGOptions*/
            std::string DGOptfileName("DGOptions.txt");
            std::string DGOptLoc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
            std::string DGOptkeyWordsDouble[2] = { "CourantNumber", "totalTime(s)" }, DGOptkeyWordsInt[5] = {"numberOfGaussPoints1D","numberOfGaussPoints2D","orderOfAccuracy", "writeInterval","totalProcess"}, DGOptkeyWordsBool[2] = { "writeLog", "loadSavedCase"}, DGOptkeyWordsStr[2] = {"ddtScheme", "runningMode"};
            double DGOptoutDB[2] = {};
            int DGOptoutInt[5] = {};
            bool DGOptoutBool[2] = {};
            std::string DGOptoutStr[2] = {};

            readDataFile(DGOptfileName, DGOptLoc, DGOptkeyWordsDouble, DGOptkeyWordsInt, DGOptkeyWordsBool, DGOptkeyWordsStr, DGOptoutDB, DGOptoutInt, DGOptoutBool, DGOptoutStr, 2, 5, 2, 2);

            systemVar::CFL = DGOptoutDB[0];
            systemVar::Ttime = DGOptoutDB[1];
            mathVar::nGauss1D = DGOptoutInt[0];
            mathVar::nGauss2D = DGOptoutInt[1];
            mathVar::orderElem = DGOptoutInt[2];
            systemVar::wrtI = DGOptoutInt[3];
            systemVar::totalProc = DGOptoutInt[4];
            systemVar::wrtLog = DGOptoutBool[0];
            systemVar::loadSavedCase = DGOptoutBool[1];

            if (DGOptoutStr[0].compare("Euler") == 0)
            {
                systemVar::ddtScheme = 1;
            }
            /*
            else if (DGOptoutStr[0].compare("TVDRK2") == 0)
            {
                systemVar::ddtScheme = 2;
            }*/
            else if (DGOptoutStr[0].compare("TVDRK3") == 0)
            {
                systemVar::ddtScheme = 3;
            }
            else {
                std::string str0("ddtScheme '"+DGOptoutStr[0]+"' is not a time discretization scheme.");
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, str0);
            }

            if (DGOptoutStr[1].compare("parallel") == 0)
            {
                systemVar::parallelMode = true;
                systemVar::readWriteMode = "p";
            }
            else if (DGOptoutStr[1].compare("sequence") == 0)
            {
                systemVar::parallelMode = false;
                systemVar::readWriteMode = "s";
                systemVar::totalProc=1;
            }
        }

        /**
         * @brief Function reads Material file
         */
        void Material()
        {
            /*Read Material*/
            std::string MatfileName("Material.txt");
            std::string MatLoc;
            std::string visModel;
            MatLoc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant";
            std::string MatkeyWordsDouble[3] = { "gammaRatio", "gasConstant", "PrandtlNumber"}, MatkeyWordsInt[1] = {}, MatkeyWordsBool[1] = {}, MatkeyWordsStr[1] = {"viscosity"};
            double MatoutDB[4] = {};
            int MatoutInt[1] = {};
            bool MatoutBool[1] = {};
            std::string MatoutStr[1] = {};

            readDataFile(MatfileName, MatLoc, MatkeyWordsDouble, MatkeyWordsInt, MatkeyWordsBool, MatkeyWordsStr, MatoutDB, MatoutInt, MatoutBool, MatoutStr, 3, 0, 0, 1);
            material::gamma = MatoutDB[0];
            material::R = MatoutDB[1];
            material::Pr = MatoutDB[2];

            material::Cp = material::R*material::gamma / (material::gamma - 1);
            material::Cv = material::Cp - material::R;

            visModel=MatoutStr[0];
            //Read viscosity settings
            if (visModel.compare("Sutherland")==0)
            {
                material::viscousityModel::sutherland=true;

                const int doubleArrSize(2);
                const int intArrSize(1);
                const int boolArrSize(1);
                const int strArrSize(1);

                //Intput so luong cac bien can doc o day
                int numDbl(2),
                        numInt(0),
                        numBool(0),
                        numStr(0);
                std::string
                        keyWordsDouble[doubleArrSize] = { "As", "Ts" },
                        keyWordsInt[intArrSize] = {},
                        keyWordsBool[boolArrSize] = {},
                        keyWordsStr[strArrSize] = {};

                double outDB[doubleArrSize] = {};
                int outInt[intArrSize] = {};
                bool outBool[boolArrSize] = {};
                std::string outStr[strArrSize] = {};

                IO::readDataFile(MatfileName, MatLoc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr,
                                 outDB, outInt, outBool, outStr,
                                 numDbl, numInt, numBool, numStr);
                material::viscosityCoeff::Sutherland::As=outDB[0];
                material::viscosityCoeff::Sutherland::Ts=outDB[1];
            }
            else if (visModel.compare("powerLaw_VHS")==0)
            {
                material::viscousityModel::power_VHS=true;

                const int doubleArrSize(4);
                const int intArrSize(1);
                const int boolArrSize(1);
                const int strArrSize(1);

                //Intput so luong cac bien can doc o day
                int numDbl(4),
                        numInt(0),
                        numBool(0),
                        numStr(0);
                std::string
                        keyWordsDouble[doubleArrSize] = { "molMass", "omega", "Tref", "dRef" },
                        keyWordsInt[intArrSize] = {},
                        keyWordsBool[boolArrSize] = {},
                        keyWordsStr[strArrSize] = {};

                double outDB[doubleArrSize] = {};
                int outInt[intArrSize] = {};
                bool outBool[boolArrSize] = {};
                std::string outStr[strArrSize] = {};

                IO::readDataFile(MatfileName, MatLoc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr,
                                 outDB, outInt, outBool, outStr,
                                 numDbl, numInt, numBool, numStr);
                material::viscosityCoeff::powerLaw_VHS::molMass=outDB[0];
                material::viscosityCoeff::powerLaw_VHS::omega=outDB[1];
                material::viscosityCoeff::powerLaw_VHS::TRef=outDB[2];
                material::viscosityCoeff::powerLaw_VHS::dRef=outDB[3];
            }
            else if (visModel.compare("constant")==0)
            {
                material::viscousityModel::constant=true;

                const int doubleArrSize(1);
                const int intArrSize(1);
                const int boolArrSize(1);
                const int strArrSize(1);

                //Intput so luong cac bien can doc o day
                int numDbl(1),
                        numInt(0),
                        numBool(0),
                        numStr(0);
                std::string
                        keyWordsDouble[doubleArrSize] = { "mu" },
                        keyWordsInt[intArrSize] = {},
                        keyWordsBool[boolArrSize] = {},
                        keyWordsStr[strArrSize] = {};

                double outDB[doubleArrSize] = {};
                int outInt[intArrSize] = {};
                bool outBool[boolArrSize] = {};
                std::string outStr[strArrSize] = {};

                IO::readDataFile(MatfileName, MatLoc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr,
                                 outDB, outInt, outBool, outStr,
                                 numDbl, numInt, numBool, numStr);
                material::viscosityCoeff::constant::mu=outDB[0];
            }
        }

        void FlowProperties()
        {
            /*Read FlowProperties*/
            std::string fileName("FlowProperties.txt");
            std::string Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
            std::string keyWordsDouble[1] = {}, keyWordsInt[1] = {}, keyWordsBool[2] = {"Viscosity", "SelfDiffusion"}, keyWordsStr[1] = {};
            double outDB[1] = {};
            int outInt[1] = {};
            bool outBool[2] = {};
            std::string outStr[1] = {};

            readDataFile(fileName, Loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr, outDB, outInt, outBool, outStr, 0, 0, 2, 0);

            flowProperties::viscous=outBool[0];
            flowProperties::massDiffusion=outBool[1];

            /*Neu set dong la inviscid, modify As va Ts bang 0*/
            if (!flowProperties::viscous)
            {
                //material::viscosityCoeff::Sutherland::As=0.0;
                //material::viscosityCoeff::Sutherland::Ts=0.0;
                /*Khong turn on massDiffusion cho dong inviscid*/
                if (flowProperties::massDiffusion)
                {
                    flowProperties::massDiffusion=false;
                    if (systemVar::currentProc==0)
                        std::cout<<"Warning! DGSolver detected that flow properties are set to inviscid and Self diffusion ON. This is unphysical and DGSolver will turn off mass diffusion.\n";
                }
            }

            else
            {
                if (flowProperties::massDiffusion)
                {
                    IO::loadSettingFiles::massDiffSettings();
                }
            }
        }

        void DGSchemes()
        {
            //Reset flags
            DGSchemes::fluxControl::LxF=false;
            DGSchemes::fluxControl::Roe=false;
            DGSchemes::fluxControl::HLL=false;
            DGSchemes::fluxControl::HLLC=false;
            DGSchemes::fluxControl::central=false;

            /*Read DGSchemes*/
            std::string fileName=("DGSchemes.txt");
            std::string Loc=(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
            std::string DGSchemeskeyWordsDouble[1] = {}, DGSchemeskeyWordsInt[1] = {}, DGSchemeskeyWordsBool[1] = {}, DGSchemeskeyWordsStr[3] = {"convectiveFlux","diffusiveTermScheme"};
            double DGSchemesoutDB[1] = {};
            int DGSchemesoutInt[1] = {};
            bool DGSchemesoutBool[1] = {};
            std::string DGSchemesoutStr[2] = {};

            readDataFile(fileName, Loc, DGSchemeskeyWordsDouble, DGSchemeskeyWordsInt, DGSchemeskeyWordsBool, DGSchemeskeyWordsStr, DGSchemesoutDB, DGSchemesoutInt, DGSchemesoutBool, DGSchemesoutStr, 0, 0, 0, 2);

            //Convective Term
            if (DGSchemesoutStr[0].compare("LxF")==0)
            {
                DGSchemes::fluxControl::LxF=true;
            }
            else if (DGSchemesoutStr[0].compare("Roe")==0)
            {
                DGSchemes::fluxControl::Roe=true;
            }
            else if (DGSchemesoutStr[0].compare("HLLE")==0)
            {
                DGSchemes::fluxControl::HLL=true;
            }
            else if (DGSchemesoutStr[0].compare("HLLC")==0) //Hien tai chua cap nhat flux nay
            {
                DGSchemes::fluxControl::HLLC=true;
            }
            else if (DGSchemesoutStr[0].compare("central")==0)
            {
                DGSchemes::fluxControl::central=true;
            }
            else {
                std::string str0("convectiveFlux type "+DGSchemesoutStr[0]+"' is not available.");
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, str0);
            }

            //Diffusive term
            if (DGSchemesoutStr[1].compare("BR1")==0)
            {
                systemVar::auxVariables=1;
            }
            else if (DGSchemesoutStr[1].compare("BR2")==0)
            {
                systemVar::auxVariables=2;
            }
            else {
                std::string str0("diffusiveTermScheme '"+DGSchemesoutStr[1]+"' is not a diffusion scheme.");
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, str0);
            }
        }

        void TBounds()
        {
            //Input do dai cua cac array o day
            const int numOfInt(1),
                    numOfDouble(2),
                    numOfBool(1),
                    numOfStr(1);

            /*Read TBounds*/
            std::string fileName("TBounds.txt");
            std::string loc;
            loc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant";
            std::string keyWordsDouble[numOfDouble] = { "upperBound", "lowerBound"}, keyWordsInt[numOfInt] = {}, keyWordsBool[numOfBool] = {}, keyWordsStr[numOfStr] = {};
            double outDB[numOfDouble] = {};
            int outInt[numOfInt] = {};
            bool outBool[numOfBool] = {};
            std::string outStr[numOfStr] = {};

            /*NOTE: trong ham readDataFile, 4 argument cuoi cung la so luong variable trong file can phai doc tuong ung voi kieu du lieu:
             * double, int, bool, string.
             * Neu kieu du lieu nao khong co data can doc, them dau '-' vao phia truoc ten bien, vi du:
             * readDataFile(..., numOfDouble, -numOfInt, -numOfBool, -numOfStr) ----> chi doc bien co kieu du lieu double
            */
            readDataFile(fileName, loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr, outDB, outInt, outBool, outStr, numOfDouble, -numOfInt, -numOfBool, -numOfStr);
            limitVal::TUp = outDB[0];
            limitVal::TDwn = outDB[1];
        }

        void massDiffSettings()
        {
            //Input do dai cua cac array o day
            const int numOfInt(1),
                    numOfDouble(1),
                    numOfBool(1),
                    numOfStr(1);

            /*Read FlowProperties*/
            std::string fileName("FlowProperties.txt");
            std::string Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
            std::string keyWordsDouble[numOfDouble] = {}, keyWordsInt[numOfInt] = {}, keyWordsBool[numOfBool] = {"SelfDiffusionAtWall"}, keyWordsStr[numOfStr] = {"Model"};
            double outDB[numOfDouble] = {};
            int outInt[numOfInt] = {};
            bool outBool[numOfBool] = {};
            std::string outStr[numOfStr] = {};

            /*NOTE: trong ham readDataFile, 4 argument cuoi cung la so luong variable trong file can phai doc tuong ung voi kieu du lieu:
             * double, int, bool, string.
             * Neu kieu du lieu nao khong co data can doc, them dau '-' vao phia truoc ten bien, vi du:
             * readDataFile(..., numOfDouble, -numOfInt, -numOfBool, -numOfStr) ----> chi doc bien co kieu du lieu double
            */
            readDataFile(fileName, Loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr, outDB, outInt, outBool, outStr, -numOfDouble, -numOfInt, numOfBool, numOfStr);

            if (outStr[0].compare("Durst")==0)
            {
                extNSF_Durst::enable=true;
                IO_Durst::readMassDiffCoefModel();
            }
            else {
                std::string str0("Model '"+outStr[0]+"' of FlowProperties>SelfDiffusionSetting is not a available.");
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, str0);
            }

            extNSF_Durst::diffusionAtWall=outBool[0];
        }
    }
}
