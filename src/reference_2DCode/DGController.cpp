#include "DGController.h"
#include "DGIOLib.h"
#include "DGMessagesLib.h"
#include "DGMeshReaderLib.h"
#include "DGMath.h"
#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "CommandCheck.h"
#include "DGProcLib.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include "dynamicVarDeclaration.h"
#include "./parallelFunctions/generalParallelFuncs.h"
#include "./parallelFunctions/parallelVariables.h"

#include "./limiters/limiterController.h"

//Non Equilibrium BCs
#include "./boundaryConditions/customBCs/nonEquilibriumBCs/nonEqmBCs_GenFuncs.h"

//Extended Navier-Stokes-Fourier model
#include "./extNSFEqns/FranzDurst/DurstModel.h"

/**
 * @brief Function excutes command.
 *
 * Functions run only at processor rank 0:
 * - exportMeshToMetis
 * - testMeshPartitionResult
 * - decomposeCase
 * - checkUnvReader
 * - checkBCsHelper
 * - reSubmit
 * - mapResults
 */
void Executer()
{
    /*Functions run only at rank 0:
        exportMeshToMetis
        testMeshPartitionResult
        decomposeCase
        checkUnvReader
        checkBCsHelper
        reSubmit
        mapResults
    */

    if (controlFlag::sequence::checkUnvReader)
    {
        if (systemVar::currentProc==0)
        {
            std::string UnvReaderLoc(systemVar::wD + "/Ultilities/MeshReader/UnvToDG");
            auxUlti::openFileEXE(UnvReaderLoc);
            controlFlag::sequence::checkUnvReader=false;
        }
    }
    else if (controlFlag::parallel::checkDGRun)
    {
        if (systemVar::runPreProcess == false)
        {
            PreProcessing();
        }
		
        Processing();
        controlFlag::parallel::checkDGRun=false;
    }
    /*
    else if (postProcessKey::checkExit(cmd))
    {
        systemVar::endKey = true;
    }
    */
    else if (controlFlag::sequence::checkUnvHelper)
    {
        if (systemVar::currentProc==0)
        {
            message::UnvReaderHelp();
        }
        controlFlag::sequence::checkUnvHelper=false;
    }
    else if (controlFlag::sequence::checkBCsHelper)
    {
        if (systemVar::currentProc==0)
        {
            message::BCsHelp();
        }
        controlFlag::sequence::checkBCsHelper=false;
    }
    else if (controlFlag::sequence::reSubmit)
    {
        if (systemVar::currentProc==0)
        {
            IO::getCase();
        }
        //PreProcessing();
        controlFlag::sequence::reSubmit=false;
    }
    else if (controlFlag::parallel::mapResults)
    {
        PreProcessing();
        meshParam::calcStiffMatrixCoeffs();
        IO::importCase::importResultsFromAnotherCase();
        systemVar::initializedOrNot = true;
        controlFlag::parallel::mapResults=false;
    }
    else if (controlFlag::sequence::exportMeshToMetis)
    {
        if (systemVar::currentProc==0)
        {
            IO::loadMesh("s");
            std::cout<<"Exporting DG mesh to Metis's format\n";
            MshExporter::exportMeshToMetis();
            controlFlag::sequence::exportMeshToMetis=false;
        }
    }
    else if (controlFlag::sequence::testMeshPartitionResult)
    {
        if (systemVar::currentProc==0)
        {
            std::string key("n");
            std::cout<<"Do you want to read mesh? <y/n>: ";
            std::cin>>key;
            if ((key.compare("y") == 0) || (key.compare("Y") == 0))
            {
                IO::loadMesh("s");
            }
            std::cout<<"Creating techplot file\n";
            MshExporter::testMeshPartitionResult();
            controlFlag::sequence::testMeshPartitionResult=false;
        }
    }
    else if (controlFlag::sequence::debug_checkElement)
    {
        int input(-1);
        std::cout << "Input element ID (ID is supplied by SALOME): ";
        std::cin >> input;
        std::cout << " \n";
        debugTool::checkElemInfor(input);
        controlFlag::sequence::debug_checkElement=false;
    }
    else if (controlFlag::sequence::decomposeCase)
    {
        if (systemVar::currentProc==0)
        {
            systemVar::runDecomposeCaseFnc=true;
            /*LOAD TOTAL PROCESSES*/
            //IO::readNumberOfCores();

            /*LOAD MESH*/
            IO::loadMesh("s");

            /*PROCESS MESH*/
            MshReader::meshProcess();

            decomposePart::decomposingMesh();

            controlFlag::sequence::decomposeCase=false;
        }
    }
    else if (controlFlag::sequence::reconstructLatestTime)
    {
        if (systemVar::currentProc==0)
        {
            /*LOAD TOTAL PROCESSES*/
            //IO::readNumberOfCores();

            /*LOAD MESH*/
            IO::loadMesh("s");

            /*PROCESS MESH*/
            MshReader::meshProcess();

            //Load time
            IO::loadTime();

            reconstructLatestTime();

            controlFlag::sequence::reconstructLatestTime=false;
        }
    }

    else if (controlFlag::sequence::checkPartitionedMesh)
    {
        if (systemVar::currentProc==0)
        {
            systemVar::runCheckParMesh=true;

            /*LOAD TOTAL PROCESSES*/
            //IO::readNumberOfCores();

            /*LOAD MESH*/
            IO::loadMesh("s");

            /*PROCESS MESH*/
            MshReader::meshProcess();

            std::cout<<"    Checking partitioned mesh...\n";
            decomposePart::fixPartitionedMesh();

            controlFlag::sequence::checkPartitionedMesh=false;
        }
    }
}

/**
 * @brief Function checks whether inputted command is available
 * @param cmd: inputted command
 */
void checkCommandLine(std::string cmd)
{
    if (preProcessKey::checkUnvReader(cmd))
    {
        controlFlag::sequence::checkUnvReader=true;
    }
    else if (processKey::checkDGRun(cmd))
    {
        controlFlag::parallel::checkDGRun=true;
    }
    else if (postProcessKey::checkExit(cmd))
    {
        systemVar::endKey = true;
    }
    else if (preProcessKey::checkUnvHelper(cmd))
    {
        controlFlag::sequence::checkUnvHelper=true;
    }
    else if (preProcessKey::checkBCsHelper(cmd))
    {
        controlFlag::sequence::checkBCsHelper=true;
    }
    else if (preProcessKey::reSubmit(cmd))
    {
        controlFlag::sequence::reSubmit=true;
    }
    else if (preProcessKey::mappResults(cmd))
    {
        controlFlag::parallel::mapResults=true;
    }
    else if (preProcessKey::exportMeshToMetis(cmd))
    {
        controlFlag::sequence::exportMeshToMetis=true;
    }
    else if (preProcessKey::testMeshPartitionResult(cmd))
    {
        controlFlag::sequence::testMeshPartitionResult=true;
    }
    else if (preProcessKey::debug::checkElement(cmd))
    {
        controlFlag::sequence::debug_checkElement=true;
    }
    else if (preProcessKey::decomposeCase(cmd))
    {
        controlFlag::sequence::decomposeCase=true;
    }
    else if (preProcessKey::reconstructLatestTime(cmd))
    {
        controlFlag::sequence::reconstructLatestTime=true;
    }
    else if (preProcessKey::checkPartitionedMesh(cmd))
    {
        controlFlag::sequence::checkPartitionedMesh=true;
    }
    else
    {
        std::cout << "	ERROR: unknow command <" << cmd << ">!!!\n";
    }
}

/**
 * @brief Function does Processing
 */
void Processing()
{
    //Setup case san sang de chay
    process::prepareCaseForRunning();

	while (process::checkRunningCond())
	{
		systemVar::iterCount++;
        systemVar::savingCout++;
        if (systemVar::currentProc==0)
        {
            std::cout << "Iteration " << systemVar::iterCount << std::endl
                      << "Total time = "<<runTime<<std::endl;
        }

        //CALCULATE TIME STEP
        process::timeDiscretization::calcGlobalTimeStep();
        process::timeDiscretization::parallel::minTimeStep();

        //SOLVE TIME MARCHING
        process::timeDiscretization::solveTimeEq();

		//COMPUTE RESIDUALS
		process::timeDiscretization::globalErrorEstimate();

        //Ham ket thuc 1 iteration
        process::finalizeIteration();
	}
}

/**
 * @brief Function does Pre-Processing.
 *
 * - Read Mesh
 * - Read Boundary Condition
 * - Check information before running (require user's confirmation)
 * - Calculate mesh's parameters
 * - Resize data arrays
 */
void PreProcessing()
{
    /*LOAD CONSTANTS*/
    IO::loadSettingFiles::loadConstants();
    //limiter::IOLimiter::readSelectedLimiters();

    std::string readWriteMode;
    if (systemVar::parallelMode)
    {
        readWriteMode="p";
    }
    else {
        readWriteMode="s";
    }

	/*LOAD MESH*/
    IO::loadMesh(readWriteMode);

	/*LOAD p T U*/
    IO::loadpTU(readWriteMode);

    /*CHECK INFROMATIONS BEFORE RUNNING*/
    auxUlti::checkInforBeforeRunning();

	/*PROCESS MESH*/
	MshReader::meshProcess();

    /*RESIZE ARRAYS*/
    auxUlti::resizeRequiredArrays();
    nonEquilibriumBCs::resizeSurfaceFields();
    extNSF_Durst::resizeArrays();

    //Synch mesh data
    for (int ivertex=0;ivertex<4;ivertex++)
    {
        parallelFuncs_Gen::sendReceiveMeshData(ivertex,1,parallelBuffer::xCoor);
        parallelFuncs_Gen::sendReceiveMeshData(ivertex,2,parallelBuffer::yCoor);
    }

	/*CALCULATE JACOBIAN, BASIS FUNCTION AND GAUSSIAN*/
    meshParam::genBasedGaussPtsVectors();
	meshParam::GaussParam();
	meshParam::basisFcParam();
	meshParam::JacobianParam();

	/*CALCULATE CELL METRICS*/
	meshParam::calcCellMetrics();
    meshParam::calcEdgeLength();

	/*CALCULATE COORDINATES DERIVATIVES*/
	meshParam::derivCoordinates();
    auxUlti::mappingEdges();

    /*Ham tim hinh chieu vuong goc cua center den BCEdge*/
    meshParam::findNormProjectionOfCenterToBCEdge();

    //Tinh khoang cach tu centroid cua cell tai BC den BC
    meshParam::calcDistanceFromCenterToBCEdge();

    //Tinh khoang cach tu centroid cua cell tai BC den diem Gauss tren BC
    math::geometricOp::findGaussPtsOnWallParameters();

	systemVar::runPreProcess = true;
}

void PostProcessing()
{
	
}
