#include "DGMessagesLib.h"
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <mpi.h>

//include chrono and ctime to get time data
#include <chrono>
#include <ctime>

#include "./limiters/massDiffusion/massDiffusion.h"

//Durst model
#include "./extNSFEqns/FranzDurst/DurstModel.h"
#include "./extNSFEqns/FranzDurst/IO.h"

namespace message
{
	std::string headerFile()
	{
		std::string headerStr(" ");
		headerStr = R"(+------------------------------ CFD LOVES C++ ---------------------------+
|--------------------DISCONTINUOS GALERKIN METHOD SOLVER-----------------|
|                                  Author                                |
|   Nam Ha Tran.                                                         |
|   Ver 1.01                                                             |
+------------------------------------------------------------------------+
|   This program uses Discontinous Galerkin method to solve 2D problems  |
|   on structural and unstructural mesh.                                 |
+------------------------------------------------------------------------+
)";
		std::string Timestr(getTime());
		headerStr += "	Program ran at (d-m-y_h-m-s) " + Timestr + "\n";
		return headerStr;
	}

	std::string undfKeyW(std::string keyW, std::string location)
	{
        std::string str("Cannot find key word <" + keyW + "> in " + location);
		return str;
	}

	std::string opFError(std::string fileName, std::string location)
	{
		std::string str("Cannot open file <" + fileName + "> located at " + location);
		return str;
	}

	void writeLog(std::string location, std::string caseName, std::string str)
	{
        if (systemVar::currentProc==0)
        {
            std::string strTime(getTime());
            std::string logFile(location + "\\log_" + caseName + ".txt");
            std::ofstream logFlux(logFile.c_str(), std::ios_base::app | std::ios_base::out);

            //Get time
            std::string crash_time(message::getTime());

            std::cout << "ERROR: " << str << std::endl;
            if (systemVar::wrtLog==true)
            {
                logFlux << "log was created at " << crash_time << std::endl << str << std::endl;
            }
            std::cout << "DGSolver will exit after you hit return.\n";
        }
		exit(EXIT_FAILURE);
	}

	std::string getTime() //This function is referenced from internet
	{
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);

		std::ostringstream oss;
		oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
		auto str = oss.str();

		return str;
	}

	std::string undfBcType(std::string keyW, std::string fileName, std::string bcType)
	{
		std::string str("Boundary condition <" + keyW + "> of file " + fileName + " is not valid for type <" + bcType + ">");
		return str;
	}

	std::string headerpTU(std::string file)
	{
		std::string headerStr(" ");
		if (file.compare("p") == 0)
		{
			headerStr = R"(Pressure conditions (Pa)
initialValue				0
)";
		}
		else if (file.compare("T") == 0)
		{
			headerStr = R"(Temperature conditions (K)
initialValue				0
)";
		}
		else if (file.compare("U") == 0)
		{
			headerStr = R"(Velocity conditions (m/s)
initialValue			0 0 0
)";
		}
		return headerStr;
	}

	std::string SlipBcCompatibleError(int bcGrp)
	{
		std::string Str(" ");
		std::string bcGrpStr(std::to_string(bcGrp));
		Str = "Boundary condition of group " + bcGrpStr + " in file U is slip, but in file T it is not temperatureJump. This is considered as a fatal error.";
		return Str;
	}

	std::string nGaussOrderElemError()
	{
		std::string Str("orderElem should be less or equal to nGauss^2");
		return Str;
	}

	std::string BcCompatibleError(int edgeGrp)
	{
		std::string Str(" ");
		std::string bcGrpStr(std::to_string(edgeGrp));
		Str = "Boundary incompatibility error is detected, please check condition of boundary group " + bcGrpStr + " in folder /0/. This is considered as a fatal error.";
		return Str;
	}

	/*Help functions----------------------------------------------------------------*/

	void UnvReaderHelp()
	{
		std::string Str(" ");
		Str = R"(	UnvReader ultility's help:
To convert unv mesh format to DG2D readable format, do following task step by step:
	- Put mesh file into case folder.
	- Open DG2D solver, enter command <UnvToDG> to console.
	- Input mesh file name to console.
	- Edit boundary conditions, system settings and run simulation.
)";
		std::cout << Str;
	}

	void BCsHelp()
	{
		std::string Str(" ");
		Str = R"(DG Solver supports the following boundary conditions:
    +-------------------+-------------------+-------------------+
    |U                  |T                  |p                  |
    +-------------------+-------------------+-------------------+
    |1. inFlow          |1. inFlow          |1. inFlow          |
    |   value u v w     |   value T         |   value p         |
    +-------------------+-------------------+-------------------+
    |2. noSlip          |2. WallIsothermal  |2. zeroGradient    |
    |                   |   Value T         |                   |
    +-------------------+-------------------+-------------------+
    |2. noSlip          |3. WallAdiabatic   |2. zeroGradient    |
    +-------------------+-------------------+-------------------+
    |4. outFlow         |4. outFlow         |4. outFlow         |
    |   value u v w     |   value T         |   value p         |
    +-------------------+-------------------+-------------------+
    |7. symmetry        |7. symmetry        |7. symmetry        |
    +-------------------+-------------------+-------------------+
    |10. matched        |10. matched        |10. matched        |
    +-------------------+-------------------+-------------------+
)";
		std::cout << Str;
	}

    void showCaseInformations()
    {
        std::string runOrNot(" ");
        std::cout<<"\nCase settings:\n";
        std::cout<<"    - Order of accuracy: "<<mathVar::orderElem<<".\n";
        std::cout<<"    - Number of Gauss points for surface integral: "<<mathVar::nGauss1D+1<<".\n";
        std::cout<<"    - Number of Gauss points for volume integral: "<<(mathVar::nGauss2D+1)*(mathVar::nGauss2D+1)<<".\n";

        std::cout<<"\nFlow properties:\n";
        if (flowProperties::subsonic)
        {
            std::cout<<"    - Flow is subsonic, Mach number is "<<flowProperties::Mach<<".\n";
        }
        else {
            std::cout<<"    - Flow is supersonic, Mach number is "<<flowProperties::Mach<<".\n";
        }

        if (flowProperties::viscous)
        {
            std::cout<<"    - Flow is viscous. Viscosity model is ";
            if (material::viscousityModel::sutherland)
            {
                std::cout<<"Sutherland.\n";
                std::cout<<"       + Coefficients:\n"
                        <<"           As                "<<material::viscosityCoeff::Sutherland::As<<"\n"
                        <<"           Ts                "<<material::viscosityCoeff::Sutherland::Ts<<"\n"
                        <<"\n";
            }
            else if (material::viscousityModel::power_VHS)
            {
                std::cout<<"Power Law (bases on Very-Hard-Sphere model).\n";
                std::cout<<"       + Coefficients:\n"
                        <<"           molMass           "<<material::viscosityCoeff::powerLaw_VHS::molMass<<" (g/mol)\n"
                        <<"           omega             "<<material::viscosityCoeff::powerLaw_VHS::omega<<"\n"
                        <<"           TRef              "<<material::viscosityCoeff::powerLaw_VHS::TRef<<" (K)\n"
                        <<"           dRef              "<<material::viscosityCoeff::powerLaw_VHS::dRef<<" (m)\n"
                        <<"\n";
            }
            else if (material::viscousityModel::constant)
            {
                 std::cout<<"Constant Viscosity.\n";
                 std::cout<<"       + Coefficients:\n"
                        <<"           mu                "<<material::viscosityCoeff::constant::mu<<"\n";
            }
        }
        else {
            std::cout<<"    - Flow is inviscid.\n";
        }

        if (flowProperties::massDiffusion)
        {
            std::cout<<"    - Self diffusion is on:\n"
                    <<"       + Type of Extened NEF Equations: ";
            if (extNSF_Durst::enable)
            {
                std::cout<<"Durst model.\n";
                IO_Durst::showSettings();
            }
        }
        else {
            std::cout<<"    - Self diffusion is off.\n";
        }

        std::cout<<"\nBoundary conditions settings:\n";
        if (auxUlti::checkTimeVaryingBCAvailable())
        {
            std::cout<<"    - Time varying BCs are set.\n";
        }
        else {
            std::cout<<"    - All BCs are not varied with time.\n";
        }

        std::cout<<"\nScheme settings:\n";
        std::cout<<"    - Time discretization scheme: ";
        if (systemVar::ddtScheme==1)
            std::cout<<" Euler.\n";
        else if (systemVar::ddtScheme==3)
            std::cout<<" Total variation diminishing Runge-Kutta order 3.\n";


        //Show flux type
        std::cout<<"    - Flux schemes: \n";
        std::cout<<"       + Flux type of auxiliary flux: Central flux.\n";
        std::cout<<"       + Flux type of convective flux: ";
        if (DGSchemes::fluxControl::LxF)
        {
            std::cout<<"Lax-Friedrichs flux.\n";
        }
        else if (DGSchemes::fluxControl::Roe)
        {
            std::cout<<"Roe Average flux.\n";
        }
        else if (DGSchemes::fluxControl::HLL)
        {
            std::cout<<"Harten-Lax-van Leer-Einfeldt (HLLE) flux.\n";
        }
        else if (DGSchemes::fluxControl::HLLC) //Hien tai chua cap nhat flux nay
        {
            std::cout<<"Harten-Lax-van Leer-Contact (HLLC) flux.\n";
        }
        else if (DGSchemes::fluxControl::central)
        {
            std::cout<<"Central flux.\n";
        }
        std::cout<<"       + Flux type of diffusive flux: Central flux.\n";


        std::cout<<"\nLimiter settings:\n";
        std::cout << "    - Selected limiter(s): ";
        if (limitVal::PositivityPreserving)
            std::cout<<" (Positivity Preserving)";
        if (limitVal::PAdaptive)
            std::cout<<" (P Adaptive)";
        if (limitVal::massDiffusion)
            std::cout<<" (Mass Diffusion Based)";
        std::cout << "\n";

        if (limitVal::massDiffusion||limitVal::PAdaptive||limitVal::PositivityPreserving)
        {
            std::cout << "    - Limiter(s) setting: \n";
            if (limitVal::PositivityPreserving)
            {
                std::cout<<"       + Positivity Preserving:\n";
                if (limitVal::PositivityPreservingSettings::version==1)
                    std::cout<<"           Version full\n";
                else
                    std::cout<<"           Version           simplified\n";
            }
            if (limitVal::massDiffusion)
            {
                std::cout<<"       + Positivity Preserving:\n"
                        <<"           maxIter           "<<limiter::massDiffusion::maxIter<<"\n"
                        <<"           Co                "<<limiter::massDiffusion::pseudoCo<<"\n"
                        <<"           DmCoeff           "<<limiter::massDiffusion::DmCoeff<<"\n"
                        <<"\n";
            }
        }
    }

    void showWarning()
    {
        //Reversed flow
        bool showWarning(false);
        int temp;
        if (systemVar::currentProc==0)
        {
            //Check o processor 0
            if (warningFlag::reversedFlowOccur)
                showWarning=true;

            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&temp, 1, MPI_INT, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (!showWarning && temp==1)
                {
                    showWarning=true;
                }
            }
            if (showWarning)
            {
                std::cout<<"Warning! Reversed flow is occured at outlet.\n";
            }
        }
        else {
            if (warningFlag::reversedFlowOccur) temp=1;
            MPI_Send(&temp, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        warningFlag::reversedFlowOccur=false;
    }
}

void exitDG(std::string str)
{
	std::cout << "ERROR: " << str << ". This is considered as fatal error." << std::endl;
    std::cout << "DGSolver is exitting.\n";
	exit(EXIT_FAILURE);
}
