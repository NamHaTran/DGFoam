/*! \brief Main DG2D file.
 */

#include "DGIOLib.h"
#include "DGController.h"
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <mpi.h>
#include "./parallelFunctions/generalParallelFuncs.h"

using namespace std;

/**
 * @brief Main function of DG2DSolver
 * @return
 */
int main()
{
    parallelFuncs_Gen::prepareParallelCase();

	/*DISPLAY LOGO*/
    if (systemVar::currentProc==0)
    {
        IO::dispLogo();
        std::cout << "                     Welcome to DG2D solver console!!\n";
    }
    IO::getCase();

    //Phai doc so cores truoc tien
    IO::readNumberOfCores();

	while (systemVar::endKey==false)
	{
        auxUlti::getCommand();

        checkCommandLine(systemVar::cmd);

        MPI_Barrier(MPI_COMM_WORLD);
        Executer();
	}
    MPI_Finalize();
	return 0;
}
