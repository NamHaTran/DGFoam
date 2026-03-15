#include "DGLimiterLib.h"
#include <iostream>
#include <vector>
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include <tuple>
#include "DGMath.h"
#include "DGAuxUltilitiesLib.h"
#include <algorithm>
#include "DGPostProcessLib.h"
#include <math.h>
#include "DGProcLib.h"
#include "./boundaryConditions/DGBCsLib.h"
#include <mpi.h>
#include "./parallelFunctions/cellData.h"

namespace limiter
{
    //ham limiter chay sau khi hoan thanh 1 step TVDRK3
    void limiter_1Step()
	{
        if (limitVal::massDiffusion && limitVal::runningMassDiffLimiter==false)
        {
            massDiffusionLimiter::limiter();
        }

        if (mathVar::orderElem != 0)
		{
            pAdaptive::limiter();

            Pp::limiter();
        }
	}

    /**
     * @brief Function limits Rho by applying Positivity Preserving Limiter.
     */
    void limitRho_PositivityPreserving()
    {
        if (mathVar::orderElem!=0)
        {
            if (limitVal::PositivityPreserving)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {
                    theta1Arr[nelem] = limiter::Pp::calcTheta1(nelem);
                }
            }
        }
    }

	namespace Pp
	{
        void limiter()
        {
            if (limitVal::PositivityPreserving)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {
                    if (!flowProperties::massDiffusion)
                    {
                        theta1Arr[nelem] = limiter::Pp::calcTheta1(nelem);
                    }
                    theta2Arr[nelem] = limiter::Pp::calcTheta2(nelem);
                }

                if (systemVar::currentProc==0)
                {
                    int numOfTotalLimitedCell(limitVal::numOfLimitCell), recvNum;
                    for (int irank=1;irank<systemVar::totalProc;irank++) {
                        MPI_Recv(&recvNum, 1, MPI_INT, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        numOfTotalLimitedCell+=recvNum;
                    }
                    if (numOfTotalLimitedCell > 0)
                    {
                        std::cout << "Posivity preserving limiter is applied at " << numOfTotalLimitedCell << " cell(s)\n";
                    }
                }
                else {
                    MPI_Send(&limitVal::numOfLimitCell, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
                limitVal::numOfLimitCell = 0;
            }
        }

        void initialiseThetaVector()
		{
			for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
			{
				theta1Arr[nelem] = 1.0;
				theta2Arr[nelem] = 1.0;
			}
		}

        double calcTheta2(int element)
        {
            //Find theta2
            double meanRho(0.0), theta1(theta1Arr[element]), theta2(0.0), meanRhoe(0.0), minRhoe(0.0);

            std::tie(minRhoe, meanRhoe) = limiter::mathForLimiter::Pp::calcMinMeanRhoe(element, theta1);
            theta2 = limiter::mathForLimiter::Pp::calcTheta2Coeff(meanRhoe, minRhoe, meanRho);

            if ((theta2 < 1))
            {
                limitVal::numOfLimitCell++;
            }
            //Reset limit flag
            limitVal::limitFlagLocal = false;
            return theta2;
        }

        double calcTheta1(int element)
        {
            double meanRho(0.0), minRho(0.0), theta1(0.0);

            //Find theta1
            minRho = limiter::mathForLimiter::Pp::calcMinRho(element);
            meanRho = rho[element][0];

            //Compute theta1
            theta1 = limiter::mathForLimiter::Pp::calcTheta1Coeff(minRho, meanRho);
            return theta1;
        }
	}

	namespace pAdaptive
	{
        void limiter()
        {
            //p-Adaptive
            if (limitVal::PAdaptive)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {

                    for (int valType = 1; valType <= 4; valType++)
                    {
                        //p-Adaptive limiter
                        limiter::pAdaptive::limiter_1Elem(nelem, valType);
                    }
                    if (limitVal::pAdaptive::limitFlagLocal == true)
                    {
                        limitVal::pAdaptive::numOfLimitCell++;
                        limitVal::pAdaptive::limitFlagLocal = false;
                    }
                }
                if (limitVal::pAdaptive::numOfLimitCell > 0)
                {
                    std::cout << "P-adaptive limiter is applied at " << limitVal::pAdaptive::numOfLimitCell << " cell(s)\n";
                    limitVal::pAdaptive::numOfLimitCell = 0;
                }
            }
        }

        void limiter_1Elem(int element, int valType)
		{
			int elemType(auxUlti::checkType(element)),
				elemJPlus(0), elemJMinus(0),  //means j+1, j-1
				elemIPlus(0), elemIMinus(0);  //means i+1, i-1
			switch (elemType)
			{
			case 3:
			{
				std::vector<double> outputVar(mathVar::orderElem, 0.0);
				elemIPlus = meshVar::neighboringElements[element][1];
				elemIMinus = meshVar::neighboringElements[element][2];
				if (mathVar::orderElem > 1)
				{
					elemJMinus = meshVar::neighboringElements[element][0];
				}

				if ((elemIPlus >= 0) & (elemIMinus >= 0) & (elemJMinus >= 0))
				{
					switch (valType)
					{
					case 1:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rho[elemIPlus][0], rho[elemIMinus][0], rho[elemJMinus][0]);
						rho[element][1] = outputVar[0];
						rho[element][2] = outputVar[1];
					}
					break;
					case 2:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rhou[elemIPlus][0], rhou[elemIMinus][0], rhou[elemJMinus][0]);
						rhou[element][1] = outputVar[0];
						rhou[element][2] = outputVar[1];
					}
					break;
					case 3:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rhov[elemIPlus][0], rhov[elemIMinus][0], rhov[elemJMinus][0]);
						rhov[element][1] = outputVar[0];
						rhov[element][2] = outputVar[1];
					}
					break;
					case 4:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rhoE[elemIPlus][0], rhoE[elemIMinus][0], rhoE[elemJMinus][0]);
						rhoE[element][1] = outputVar[0];
						rhoE[element][2] = outputVar[1];
					}
					break;
					default:
						break;
					}
				}
			}
			break;
			case 4:
			{
				std::vector<double> outputVar(mathVar::orderElem, 0.0);
                //bool internalElem(true);
				elemIPlus = meshVar::neighboringElements[element][1];
				elemIMinus = meshVar::neighboringElements[element][3];
				if (mathVar::orderElem > 1)
				{
					elemJPlus = meshVar::neighboringElements[element][2];
					elemJMinus = meshVar::neighboringElements[element][0];
				}

				if ((elemIPlus >= 0) & (elemIMinus >= 0) & (elemJPlus >= 0) & (elemJMinus >= 0))
				{
					switch (valType)
					{
					case 1:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rho[elemIPlus][0], rho[elemIMinus][0], rho[elemJPlus][0], rho[elemJMinus][0]);
						rho[element][1] = outputVar[0];
						rho[element][2] = outputVar[1];
					}
					break;
					case 2:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rhou[elemIPlus][0], rhou[elemIMinus][0], rhou[elemJPlus][0], rhou[elemJMinus][0]);
						rhou[element][1] = outputVar[0];
						rhou[element][2] = outputVar[1];
					}
					break;
					case 3:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rhov[elemIPlus][0], rhov[elemIMinus][0], rhov[elemJPlus][0], rhov[elemJMinus][0]);
						rhov[element][1] = outputVar[0];
						rhov[element][2] = outputVar[1];
					}
					break;
					case 4:
					{
						outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rhoE[elemIPlus][0], rhoE[elemIMinus][0], rhoE[elemJPlus][0], rhoE[elemJMinus][0]);
						rhoE[element][1] = outputVar[0];
						rhoE[element][2] = outputVar[1];
					}
					break;
					default:
						break;
					}
				}
			}
			break;
			default:
				break;
			}
		}

		std::vector<double> pAdaptiveChildFunction_Quad(int element, int valType, double IPlus, double IMinus, double JPlus, double JMinus)
		{
			std::vector<double>elementConsValOfOrder(mathVar::orderElem + 1, 0.0),
				inputArgument(3, 0.0), output(mathVar::orderElem, 0.0);
			elementConsValOfOrder = auxUlti::getElementConserValuesOfOrder(element, valType);

			double UBC(math::pointValueNoLimiter(element, 1.0, 0.0, valType) - elementConsValOfOrder[0]), 
				UAD(math::pointValueNoLimiter(element, -1.0, 0.0, valType) - elementConsValOfOrder[0]),
				UCD(math::pointValueNoLimiter(element, 0.0, 1.0, valType) - elementConsValOfOrder[0]),
				UAB(math::pointValueNoLimiter(element, 0.0, -1.0, valType) - elementConsValOfOrder[0]),
				UBCMod(0.0), UADMod(0.0), UCDMod(0.0), UABMod(0.0),
				UBC_check(0.0), UAD_check(0.0), UCD_check(0.0), UAB_check(0.0),
				M(limiter::mathForLimiter::calM(element, valType)), Lxy(meshVar::localCellSize[element]);

			inputArgument[0] = UBC;
			inputArgument[1] = IPlus - elementConsValOfOrder[0];
			inputArgument[2] = elementConsValOfOrder[0] - IMinus;
			UBCMod = limiter::mathForLimiter::minmod(inputArgument);
			UBC_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1); //*fabs(elementConsValOfOrder[0])

			inputArgument[0] = -UAD;
			UADMod = -limiter::mathForLimiter::minmod(inputArgument);
			inputArgument[0] = UAD;
			UAD_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);  //*fabs(elementConsValOfOrder[0])
            if (fabs(UBC - UBC_check)>1e-10 || fabs(UAD - UAD_check)>1e-10)
			{
				//std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
				if (limitVal::pAdaptive::limitFlagLocal == false)
				{
					limitVal::pAdaptive::limitFlagLocal = true;
				}
				output[0] = 0.5*(UBCMod - UADMod);
			}
			else
			{
				output[0] = elementConsValOfOrder[1];
			}

			if (mathVar::orderElem > 1)
			{
				inputArgument[0] = UCD;
				inputArgument[1] = JPlus - elementConsValOfOrder[0];
				inputArgument[2] = elementConsValOfOrder[0] - JMinus;
				UCDMod = limiter::mathForLimiter::minmod(inputArgument);
				UCD_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);

				inputArgument[0] = -UAB;
				UABMod = -limiter::mathForLimiter::minmod(inputArgument);
				inputArgument[0] = UAB;
				UAB_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);
                if (fabs(UCD - UCD_check)>1e-10 || fabs(UAB - UAB_check)>1e-10)
				{
					//std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
					output[1] = 0.5*(UCDMod - UABMod);
				}
				else
				{
					output[1] = elementConsValOfOrder[2];
				}
			}
			return output;
		}

		std::vector<double> pAdaptiveChildFunction_Tri(int element, int valType, double IPlus, double IMinus, double JMinus)
		{
			std::vector<double>elementConsValOfOrder(mathVar::orderElem + 1, 0.0),
				inputArgument(3, 0.0), output(mathVar::orderElem, 0.0);
			elementConsValOfOrder = auxUlti::getElementConserValuesOfOrder(element, valType);

			double UBC(math::pointValueNoLimiter(element, 1.0, 0.0, valType) - elementConsValOfOrder[0]),
				UAD(math::pointValueNoLimiter(element, -1.0, 0.0, valType) - elementConsValOfOrder[0]),
				UAB(math::pointValueNoLimiter(element, 0.0, -1.0, valType) - elementConsValOfOrder[0]),
				UBCMod(0.0), UADMod(0.0), UABMod(0.0),
				UBC_check(0.0), UAD_check(0.0), UAB_check(0.0),
				M(limiter::mathForLimiter::calM(element, valType)), Lxy(meshVar::localCellSize[element]);

			inputArgument[0] = -UAB;
			inputArgument[1] = elementConsValOfOrder[0] - JMinus;
			inputArgument[2] = elementConsValOfOrder[0] - JMinus;
			UABMod = -limiter::mathForLimiter::minmod(inputArgument);
			//inputArgument[0] = UAB;
			UAB_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);
            if (fabs(-UAB - UAB_check)>1e-10)
			{
				//std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
				output[0] = -UABMod;
			}
			else
			{
				output[0] = elementConsValOfOrder[1];
			}

			if (mathVar::orderElem > 1)
			{
				inputArgument[0] = UBC;
				inputArgument[1] = IPlus - elementConsValOfOrder[0];
				inputArgument[2] = elementConsValOfOrder[0] - IMinus;
				UBCMod = limiter::mathForLimiter::minmod(inputArgument);
				UBC_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1); //*fabs(elementConsValOfOrder[0])

				inputArgument[0] = -UAD;
				UADMod = -limiter::mathForLimiter::minmod(inputArgument);
				//inputArgument[0] = UAD;
				UAD_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);  //*fabs(elementConsValOfOrder[0])
                if ((fabs(UBC - UBC_check)>1e-10 || fabs(-UAD - UAD_check)>1e-10))
				{
					//std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
					if (limitVal::pAdaptive::limitFlagLocal == false)
					{
						limitVal::pAdaptive::limitFlagLocal = true;
					}
					output[1] = (UBCMod - UADMod) / 2.0;
				}
				else
				{
					output[1] = elementConsValOfOrder[2];
				}
			}
			return output;
		}
	}

    namespace massDiffusionLimiter {
        void limiter()
        {
            limitVal::runningMassDiffLimiter=true;
            bool limit(false);
            int iter(0);

            //save current time step
            double dt_temp(dt);

            //modify dt
            dt*=0.05;

            limit = massDiffusionLimiter::checkTroubleCells();
            if (limit)
            {
                do {
                    std::cout <<"   -Iter "<<iter+1<<": ";
                    //Giai pt NSF mo rong bang TVDRK3
                    rho0 = rho;

                    for (int iRKOrder = 1; iRKOrder <= 3; iRKOrder++)
                    {
                        massDiffusionLimiter::mathForMassDiffLimiter::limiter_1step(iRKOrder);
                        massDiffusionLimiter::mathForMassDiffLimiter::updateVariables();
                        limiter::limiter_1Step();
                    }
                    iter++;
                }
                while (massDiffusionLimiter::checkRunning()&&iter<=100);
            }

            //tra runningMassDiffLimiter ve false
            limitVal::runningMassDiffLimiter=false;
            //return dt to original value
            dt=dt_temp;
        }

        bool checkTroubleCells()
        {
            bool trouble(false);
            int counter(0);
            for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
            {
                if (troubledCellDetection::checkTroubleCell_massDiff(nelem))
                {
                    limitVal::troubleCellsMarker[nelem]=true;
                    trouble=true;
                    counter++;
                }
                else
                    limitVal::troubleCellsMarker[nelem]=false;
            }

            if (counter!=0)
            {
                std::cout << "massDiffusion limiter is applied at " << counter << " cell(s)\n";
            }
            return trouble;
        }

        bool checkRunning()
        {
            bool run(false);
            int numTroubleCells(0);
            for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
            {
                if (mathForMassDiffLimiter::calcMinRhoResidual(nelem)>1e-5)
                {
                    run=true;
                    limitVal::troubleCellsMarker[nelem]=true;
                    numTroubleCells++;
                }
                else
                {
                    limitVal::troubleCellsMarker[nelem]=false;
                }
            }
            std::cout <<"Mass diffusion limiter is not converged at "<<numTroubleCells<<" cell.\n";
            return run;
        }

        namespace mathForMassDiffLimiter {
            void limiter_1step(int RKOrder)
            {
                for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
                {
                    if (limitVal::troubleCellsMarker[nelem])
                    {
                        massDiffusionLimiter::mathForMassDiffLimiter::calcVolumeGaussRho(nelem);

                        massDiffusionLimiter::mathForMassDiffLimiter::calcRhoAtElementEdge(nelem);

                        massDiffusionLimiter::mathForMassDiffLimiter::solveDivRho(nelem);

                        massDiffusionLimiter::mathForMassDiffLimiter::calcFinvFvisAtInterface(nelem);

                        massDiffusionLimiter::mathForMassDiffLimiter::solveMassEquation(RKOrder,nelem);
                    }
                }
            }

            void calcVolumeGaussRho(int element)
            {
                double a(0.0), b(0.0);
                std::vector<double>UVol(4,0.0);
                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss+1));
                        std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                        volumeFields::rhoVolGauss[element][nanb] = rho[element][0];
                        for (int iorder = 1; iorder <= mathVar::orderElem; iorder++)
                        {
                            volumeFields::rhoVolGauss[element][nanb] += rho[element][iorder]* mathVar::B[iorder];
                        }
                    }
                }
            }

            void calcRhoAtElementEdge(int element)
            {
                int bcGrp(0), elemType(auxUlti::checkType(element)), iedge(0), loc(0);
                double aSurf(0.0), bSurf(0.0);

                for (int i = 0; i < elemType; i++)
                {
                	//get id of edge
                	iedge=meshVar::inedel[element][i];
                    bcGrp = auxUlti::getBCType(iedge);
                    if (bcGrp == 0)
                    {
                        for (int nG = 0; nG <= mathVar::nGauss; nG++)
                        {
                            std::tie(aSurf, bSurf) = auxUlti::getGaussSurfCoor(iedge, element, nG);
                            if (auxUlti::checkMaster(element,iedge)) loc=nG;
                            else loc=nG + mathVar::nGauss + 1;

                            surfaceFields::rho[iedge][loc] = math::pointValue(element, aSurf, bSurf, 1, 2);
                        }
                    }
                }
            }

            void solveDivRho(int element)
	        {
	            std::vector<double> rhoRHSTermOxDir(mathVar::orderElem + 1, 0.0),
	                rhoRHSTermOyDir(mathVar::orderElem + 1, 0.0);

	            //2) Calculate Right hand side terms
                process::auxEq::massDiffusion::BR1::CalcRHSTerm(element, rhoRHSTermOxDir, rhoRHSTermOyDir);

                //3) Solve for div(rho)
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    //Ox direction
                    //BR1Vars::dRhoX[element][iorder] = rhoRHSTermOxDir[iorder] / stiffMatrixCoeffs[element][iorder];

                    //Oy direction
                    //BR1Vars::dRhoY[element][iorder] = rhoRHSTermOyDir[iorder] / stiffMatrixCoeffs[element][iorder];
                }
	        }

            void calcFinvFvisAtInterface(int element)
            {
                /* gia tri cua T, rhou, rhov volume va surface deu lay gia tri cua vong lap truoc nen
                 * khong can tinh lai
                */

                int masterCell(-1), slaveCell(-1), bcGrp(0), elemType(auxUlti::checkType(element)), iedge(0);
                double uMaster(0.0), vMaster(0.0), TMaster(0.0), uSlave(0.0), vSlave(0.0), TSlave(0.0),
                    uMagM(0.0), uMagP(0.0), aM(0.0), aP(0.0), muMaster(0.0), muSlave(0.0),
                    rhoMaster(0.0), dRhoXMaster(0.0), dRhoYMaster(0.0),
                    rhoSlave(0.0), dRhoXSlave(0.0), dRhoYSlave(0.0),
                    rhouMaster(0.0), rhovMaster(0.0), rhouSlave(0.0), rhovSlave(0.0);
                std::vector<double> CArray(mathVar::nGauss + 1, 0.0);

                for (int i = 0; i < elemType; i++)
                {
                    iedge=meshVar::inedel[element][i];
                    bcGrp = auxUlti::getBCType(iedge);
                    if (bcGrp == 0)
                    {
                        std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                        for (int nG = 0; nG <= mathVar::nGauss; nG++)
                        {
                            std::tie(TMaster,TSlave)=auxUlti::getTAtInterfaces(iedge,masterCell,nG);
                            muMaster=math::CalcVisCoef(TMaster);
                            muSlave=math::CalcVisCoef(TSlave);

                            std::tie(rhoMaster, rhoSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 1);
                            std::tie(rhouMaster, rhouSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 2);
                            std::tie(rhovMaster, rhovSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 3);
                            //std::tie(dRhoXMaster, dRhoXSlave) = math::internalSurfaceDerivRhoValue(iedge, masterCell, nG, 1);
                            //std::tie(dRhoYMaster, dRhoYSlave) = math::internalSurfaceDerivRhoValue(iedge, masterCell, nG, 2);
                            dRhoXMaster*=muMaster;
                            dRhoYMaster*=muMaster;
                            dRhoXSlave*=muSlave;
                            dRhoYSlave*=muSlave;

                            uMaster = (rhouMaster / rhoMaster);
                            uSlave = (rhouSlave / rhoSlave);

                            vMaster = (rhovMaster / rhoMaster);
                            vSlave = (rhovSlave / rhoSlave);

                            /*INVISCID TERMS*/
                            /*calculate velocity magnitude*/
                            uMagP = sqrt(pow(uMaster, 2) + pow(vMaster, 2));
                            uMagM = sqrt(pow(uSlave, 2) + pow(vSlave, 2));

                            /*calculate speed of sound*/
                            aP = math::CalcSpeedOfSound(TMaster);
                            aM = math::CalcSpeedOfSound(TSlave);

                            /*calculate constant for Lax-Friederich flux*/
                            CArray[nG] = math::numericalFluxes::constantC(uMagP, uMagM, aP, aM);

                            /*calculate inviscid terms*/
                            //surfaceFields::invis_rhoX[iedge][nG] = rhoMaster * uMaster;
                            //surfaceFields::invis_rhoX[iedge][nG + mathVar::nGauss + 1] = rhoSlave * uSlave;
                            //surfaceFields::invis_rhoY[iedge][nG] = rhoMaster * vMaster;
                            //surfaceFields::invis_rhoY[iedge][nG + mathVar::nGauss + 1] = rhoSlave * vSlave;

                            /*calculate viscous terms*/
                            //surfaceFields::Vis_rhoX[iedge][nG]=-material::massDiffusion::DmCoeff*dRhoXMaster/rhoMaster;
                            //surfaceFields::Vis_rhoX[iedge][nG + mathVar::nGauss + 1]=-material::massDiffusion::DmCoeff*dRhoYMaster/rhoMaster;
                            //surfaceFields::Vis_rhoY[iedge][nG]=-material::massDiffusion::DmCoeff*dRhoYMaster/rhoMaster;
                            //surfaceFields::Vis_rhoY[iedge][nG + mathVar::nGauss + 1]=-material::massDiffusion::DmCoeff*dRhoYMaster/rhoMaster;
                        }
                        LxFConst[iedge] = *std::max_element(CArray.begin(), CArray.end());
                        //DiffusiveFluxConst[iedge] = *std::max_element(BetaArray.begin(), BetaArray.end());
                    }
                }
            }

            std::tuple<double, double, double, double> calcFinvFvisInVolume(int element, int na, int nb)
            {
                int nanb(calcArrId(na,nb,mathVar::nGauss+1));
                double
                    rhoVal(volumeFields::rhoVolGauss[element][nanb]),
                    rhouVal(volumeFields::rhouVolGauss[element][nanb]),
                    rhovVal(volumeFields::rhovVolGauss[element][nanb]),
                    uVal(0.0),
                    vVal(0.0),
                    a(0.0),b(0.0),muVal(muVal=math::CalcVisCoef(volumeFields::T[element][nanb])),
                    dRhoX(0.0), dRhoY(0.0),
                    invTermX(0.0), visTermX(0.0), invTermY(0.0), visTermY(0.0);

                std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                uVal = rhouVal / rhoVal;
                vVal = rhovVal / rhoVal;
                dRhoX=math::pointAuxValue(element,a,b,1,1)*muVal;
                dRhoY=math::pointAuxValue(element,a,b,1,2)*muVal;

                /*1. Ox direction*/
                invTermX=rhoVal*uVal;
                visTermX=-material::massDiffusion::DmCoeff*dRhoX/rhoVal;

                /*2. Oy direction*/
                invTermX=rhoVal*vVal;
                visTermX=-material::massDiffusion::DmCoeff*dRhoY/rhoVal;
                return std::make_tuple(invTermX,visTermX,invTermY,visTermY);
            }

            void calcVolumeIntegralTerms(int element, std::vector<double>&VolIntTerm1)
            {
                std::vector<std::vector<double>>
                    GsVolX1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
                    GsVolY1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
                double invTermX(0.0), visTermX(0.0), invTermY(0.0), visTermY(0.0);

                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        std::tie(invTermX,visTermX,invTermY,visTermY)=mathForMassDiffLimiter::calcFinvFvisInVolume(element,na,nb);
                        /*A INVISCID TERMS*/
                        /*A1. Inviscid term on Ox direction*/
                        GsVolX1[na][nb] = invTermX;
                        /*A2. Inviscid term on Oy direction*/
                        GsVolY1[na][nb] = invTermY;

                        /*B VISCOUS TERMS*/
                        /*B1. Viscous term on Ox direction*/
                        GsVolX1[na][nb] += visTermX;
                        /*B2. Viscous term on Oy direction*/
                        GsVolY1[na][nb] += visTermY;
                    }
                }

                for (int order = 1; order <= mathVar::orderElem; order++)
                {
                    /*CALCULATE INTEGRALS*/
                    VolIntTerm1[order] = process::volumeInte(element, GsVolX1, order, 1);
                    VolIntTerm1[order] += process::volumeInte(element, GsVolY1, order, 2);
                }
                VolIntTerm1[0] = 0;
            }

            void calcSurfaceIntegralTerms(int element, std::vector<double>&SurfIntTerm1)
            {
                int elemType(auxUlti::checkType(element)), edgeName(0);
                int faceBcType(0);

                std::vector<double> Flux1Temp(mathVar::nGauss + 1, 0.0);

                std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));

                std::vector<std::vector<double>>Flux1(mathVar::nGauss + 1, std::vector<double>(4, 0.0));

                //std::vector<double> SurInt(4, 0.0);

                for (int nface = 0; nface < elemType; nface++)
                {
                    edgeName = meshVar::inedel[element][nface];
                    faceBcType = auxUlti::getBCType(edgeName);

                    if (faceBcType == 0)  //internal edge
                    {
                        for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
                        {
                            double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2));
                            double
                                termX1P(0.0), termX1M(0.0),  //(rho*u)					or 0
                                termY1P(0.0), termY1M(0.0);  //(rho*v)					or 0
                            double rhoPlus(0.0), rhoMinus(0.0);

                            /*INVISCID TERM*/
                            //Get value
                            //std::tie(termX1P, termX1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 1, 1);
                            //std::tie(termY1P, termY1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 2, 1);
                            //std::tie(rhoPlus, rhoMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1);

                            /*VISCOUS TERM*/
                            //Get value
                            //std::tie(termX1P, termX1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 1, 1);
                            //std::tie(termY1P, termY1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 2, 1);

                            /*Calculate fluxes*/
                            //Flux1[nGauss][nface] = math::numericalFluxes::advectiveFlux(termX1P, termX1M, rhoPlus, rhoMinus, LxFConst[edgeName], nx, true) + math::numericalFluxes::advectiveFlux(termY1P, termY1M, rhoPlus, rhoMinus, LxFConst[edgeName], ny, true) +
                                    //math::numericalFluxes::diffusiveFlux(termX1M, termX1P, rhoPlus, rhoMinus, DiffusiveFluxConst[edgeName], nx) + math::numericalFluxes::diffusiveFlux(termY1M, termY1P, rhoPlus, rhoMinus, DiffusiveFluxConst[edgeName], ny);
                        }
                    }
                    else  //boundary edge
                    {
                        for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
                        {
                            //Fluxes = NSFEqBCsImplement(element, edgeName, nGauss);
                            Flux1[nGauss][nface] = Fluxes[0][0] + Fluxes[0][1];
                        }
                    }
                }

                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    for (int nface = 0; nface < elemType; nface++)
                    {
                        edgeName = meshVar::inedel[element][nface];  //A BIG BUG!!!!!!
                        for (int nG = 0; nG <= mathVar::nGauss; nG++)
                        {
                            Flux1Temp[nG] = Flux1[nG][nface];
                        }
                        SurfIntTerm1[order] += process::surfaceInte(element, edgeName, Flux1Temp, order);
                    }
                }
            }

            void solveMassEquation(int RKOrder, int element)
            {
                std::vector<double> rhoError(meshVar::nelem2D, 1.0);

                std::vector<double>
                    RHSTerm1(mathVar::orderElem + 1, 0.0),
                    ddtRhoVector(mathVar::orderElem + 1, 0.0),
                    rhoVectorN(mathVar::orderElem + 1, 0.0),
                    UnVector(mathVar::orderElem + 1, 0.0);

                //2) Calculate Right hand side terms
                std::vector<double>
                    VolIntTerm1(mathVar::orderElem + 1, 0.0),
                    SurfIntTerm1(mathVar::orderElem + 1, 0.0);

                /*Volume integral term===========================================================================*/
                mathForMassDiffLimiter::calcVolumeIntegralTerms(element, VolIntTerm1);

                /*Surface integral term===========================================================================*/
                mathForMassDiffLimiter::calcSurfaceIntegralTerms(element, SurfIntTerm1);

                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    RHSTerm1[order] = VolIntTerm1[order] - SurfIntTerm1[order];
                }

                //3) Solve for time derivartives of conservative variables
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    ddtRhoVector[iorder] = RHSTerm1[iorder] / stiffMatrixCoeffs[element][iorder];
                }

                //4) Solve time marching
                //rho
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    UnVector[order] = rho[element][order];
                }
                rhoVectorN = process::NSFEq::solveTimeMarching(element, ddtRhoVector, UnVector, RKOrder, 1);

                //5) Save results to conservative variables arrays
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    rhoN[element][order] = rhoVectorN[order];
                }
            }

            void calcVolIntegralTerm(int element,std::vector<double>&VolIntTerm1)
            {
                std::vector<std::vector<double>> ViscousTerms(4, std::vector<double>(2, 0.0));
                std::vector<std::vector<double>> InviscidTerms(4, std::vector<double>(2, 0.0));
                //std::vector<double> VolInt(4, 0.0);

                std::vector<std::vector<double>>
                    GsVolX1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
                    GsVolY1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

                double rhoVal,rhouVal,rhovVal,uVal(0.0),vVal(0.0),a(0.0),b(0.0),muVal(0.0),dRhoX(0.0),dRhoY(0.0);

                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss+1));

                        //std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                        /*A INVISCID TERMS*/
                        rhoVal=volumeFields::rhoVolGauss[element][nanb];
                        rhouVal=volumeFields::rhouVolGauss[element][nanb];
                        rhovVal=volumeFields::rhovVolGauss[element][nanb];

                        std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                        uVal = rhouVal / rhoVal;
                        vVal = rhovVal / rhoVal;

                        /*A1. Inviscid term on Ox direction*/
                        GsVolX1[na][nb] = rhoVal*uVal;
                        /*A2. Inviscid term on Oy direction*/
                        GsVolY1[na][nb] = rhoVal*vVal;

                        /*B VISCOUS TERMS*/
                        muVal=math::CalcVisCoef(volumeFields::T[element][nanb]);
                        dRhoX = math::pointAuxValue(element, a, b, 1, 1)*muVal;
                        dRhoY = math::pointAuxValue(element, a, b, 1, 2)*muVal;

                        /*B1. Viscous term on Ox direction*/
                        GsVolX1[na][nb] += -material::massDiffusion::DmCoeff*dRhoX/rhoVal;
                        /*B2. Viscous term on Oy direction*/
                        GsVolY1[na][nb] += -material::massDiffusion::DmCoeff*dRhoY/rhoVal;
                    }
                }

                for (int order = 1; order <= mathVar::orderElem; order++)
                {
                    /*CALCULATE INTEGRALS*/
                    VolIntTerm1[order] = process::volumeInte(element, GsVolX1, order, 1);
                    VolIntTerm1[order] += process::volumeInte(element, GsVolY1, order, 2);
                }
                VolIntTerm1[0] = 0;
            }

            void updateVariables()
            {
                for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
                {
                    if (limitVal::troubleCellsMarker[nelement])
                    {
                        for (int order = 0; order <= mathVar::orderElem; order++)
                        {
                            rho[nelement][order] = rhoN[nelement][order];
                            rhou[nelement][order] = rhouN[nelement][order];
                            rhov[nelement][order] = rhovN[nelement][order];
                            rhoE[nelement][order] = rhoEN[nelement][order];
                        }
                    }
                }
            }

            double calcMinRhoResidual(int element)
            {
                std::vector<double> vectorRho((mathVar::nGauss+1)*(mathVar::nGauss+5), 1.0),
                        vectorRho0((mathVar::nGauss+1)*(mathVar::nGauss+5), 1.0);
                double aG(0.0), bG(0.0), minRho(0.0), minRho0(0.0);
                int counter(0);

                //Compute rho at all internal Gauss point

                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);
                        vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                        vectorRho0[counter]=(mathForMassDiffLimiter::pointRho0(element, aG, bG));
                        counter++;
                    }
                }

                //Compute rho at edge DA
                aG = -1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(mathForMassDiffLimiter::pointRho0(element, aG, bG));
                    counter++;
                }
                //Compute rho at edge BC
                aG = 1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(mathForMassDiffLimiter::pointRho0(element, aG, bG));
                    counter++;
                }
                //Compute rho at edge AB
                bG = -1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(mathForMassDiffLimiter::pointRho0(element, aG, bG));
                    counter++;
                }
                //Compute rho at edge CD
                bG = 1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(mathForMassDiffLimiter::pointRho0(element, aG, bG));
                    counter++;
                }
                minRho = *std::min_element(vectorRho.begin(), vectorRho.end());  //find min value of vector
                minRho0 = *std::min_element(vectorRho0.begin(), vectorRho0.end());
                return (fabs(minRho-minRho0)/minRho);
            }

            double pointRho0(int element, double a, double b)
            {
                double out(0.0);
                std::vector<double> Value(mathVar::orderElem + 1, 0.0);
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Value[iorder] = rho0[element][iorder];
                }

                math::basisFc(a, b, auxUlti::checkType(element));
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    out += Value[order] * mathVar::B[order];
                }

                return out;
            }
        }
    }

	//MATHEMATIC FUNCTIONS FOR LIMITER
	namespace mathForLimiter
	{
        namespace Pp {
            double calcMinRho(int element)
            {
                std::vector<double> vectorRho((mathVar::nGauss+1)*(mathVar::nGauss+5), 1.0);
                double aG(0.0), bG(0.0), min(0.0);
                int counter(0);

                //Compute rho at all internal Gauss point

                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);
                        vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                        counter++;
                    }
                }

                //Compute rho at edge DA
                aG = -1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    counter++;
                }
                //Compute rho at edge BC
                aG = 1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    counter++;
                }
                //Compute rho at edge AB
                bG = -1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    counter++;
                }
                //Compute rho at edge CD
                bG = 1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    counter++;
                }
                min = *std::min_element(vectorRho.begin(), vectorRho.end());  //find min value of vector
                return min;
            }

            //Function calculates modified value of Rho at abitrary point (for calculating theta2)
            double calcRhoModified(int element, double a, double b, double theta1)
            {
                double rhoMod(0.0);
                std::vector<double> Value(mathVar::orderElem + 1, 0.0);

                Value = auxUlti::getElementConserValuesOfOrder(element, 1);
                math::basisFc(a, b, auxUlti::checkType(element));
                for (int order = 1; order <= mathVar::orderElem; order++)
                {
                    rhoMod += Value[order] * mathVar::B[order] * theta1;
                }
                rhoMod += Value[0];
                return rhoMod;
            }

            std::tuple<double, double> calcMinMeanRhoe(int element, double theta1)
            {
                /* Neu co mass diffusion, phai modify lai rhou va rhov thanh rhou_m va rhov_m.
                 * Vi khi giai T explicit, ham CalcTFromConsvVar_massDiff tinh gia tri mass
                 * diffusion flux Jd = (1/Sc)*mu*div(rho)/rho voi Sc la so Schmidt (gia tri DmCoeff
                 * la 1/Sc) va mu tinh tu T cua step truoc (T_old), nen de limit dung can phai tinh
                 * Jd o ham limiter giong nhu o ham CalcTFromConsvVar_massDiff.
                */

                std::vector<double> vectorRhoe((mathVar::nGauss+1)*(mathVar::nGauss+5), 1.0);
                double aG(0.0), bG(0.0), min(0.0), rhoMod(0.0), rhouOrg(0.0), rhovOrg(0.0), rhoEOrg(0.0), sum(0), mean(0);
                int counter(0), edgeId;

                //Compute rho at all internal Gauss point
                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            int nanb(calcArrId(na,nb,mathVar::nGauss+1));
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(volumeFields::T[element][nanb]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod,rhouOrg,rhovOrg,rhoEOrg);
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                }

                if (auxUlti::checkType(element)==4) // cell tu giac
                {
                    //Tinh rhoe tren edge
                    int edgeAB(0), edgeBC(1), edgeCD(2), edgeDA(3);

                    //Compute rho at edge DA
                    aG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeDA);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=(limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                    //Compute rho at edge BC
                    aG = 1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeBC);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=(limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                    //Compute rho at edge AB
                    bG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        aG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeAB);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=(limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                    //Compute rho at edge CD

                    bG = 1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        aG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeCD);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=(limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                }
                else if (auxUlti::checkType(element)==3) // cell tam giac
                {
                    //Tinh rhoe tren edge
                    int edgeAB(0), edgeBC(1), edgeCA(2);

                    //Compute rho at edge CA
                    aG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeCA);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=(limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                    //Compute rho at edge BC
                    aG = 1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeBC);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=(limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                    //Compute rho at edge AB
                    bG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        aG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::mathForLimiter::Pp::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeAB);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        vectorRhoe[counter]=(limiter::mathForLimiter::Pp::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        sum+=vectorRhoe[counter];
                        counter++;
                    }
                }

                min = *std::min_element(vectorRhoe.begin(), vectorRhoe.end());  //find min value of vector
                mean = sum/counter;
                return std::make_tuple(min, mean);
            }

            double calRhoeFromConserVars(double rho, double rhou, double rhov, double rhoE)
            {
                double rhoe(0.0);
                rhoe = rhoE - 0.5*(rhou*rhou + rhov * rhov) / rho;
                return rhoe;
            }

            double calcTheta1Coeff(double minRho, double meanRho)
            {
                double temp1(0.0), theta(0.0);
                std::vector<double> vectorOmega(2, 0.0);
                vectorOmega[0] = systemVar::epsilon;
                vectorOmega[1] = meanRho;
                double omega(*std::min_element(vectorOmega.begin(), vectorOmega.end()));  //find min value of vector

                //Compare meanRho and minRho, if they are equal then theta = 1.0 (means no limiter applied)
                if (fabs(meanRho - minRho/meanRho) <= 1e-8)
                {
                    theta = 1.0;
                }
                else
                {
                    temp1 = (meanRho - omega) / (meanRho - minRho);
                    if (temp1 < 1.0)
                    {
                        theta = temp1;
                    }
                    else
                    {
                        theta = 1.0;
                    }
                }
                return theta;
            }

            double calcTheta2Coeff(double meanRhoe, double minRhoe, double meanRho)
            {
                double temp1(0.0), theta(0.0);
                std::vector<double> vectorOmega(3, 0.0);
                vectorOmega[0] = systemVar::epsilon;
                vectorOmega[1] = meanRho;
                vectorOmega[2] = meanRhoe;
                double omega(*std::min_element(vectorOmega.begin(), vectorOmega.end()));  //find min value of vector

                //Compare meanRhoe and minRhoe, if they are equal then theta = 1.0 (means no limiter applied)
                if (fabs(meanRhoe - minRhoe/meanRhoe) <= 1e-8)
                {
                    theta = 1.0;
                }
                else
                {
                    temp1 = (meanRhoe - omega) / (meanRhoe - minRhoe);
                    if (temp1 < 1.0)
                    {
                        theta = temp1;
                    }
                    else
                    {
                        theta = 1.0;
                    }
                }

                return theta;
            }
        }

		double minmod(std::vector<double> vectorArguments)
		{
			int numberArgument(vectorArguments.size());
			std::vector<double>absVectorArguments(numberArgument, 0.0);
			double minmodVal(0.0);
			bool isSameSign(true);
			for (int i = 0; i < numberArgument-1; i++)
			{
				if (limiter::mathForLimiter::getSignOfDouble(vectorArguments[i]) != limiter::mathForLimiter::getSignOfDouble(vectorArguments[i+1]))
				{
					isSameSign = false;
					break;
				}
				absVectorArguments[i] = fabs(vectorArguments[i]);
			}
			absVectorArguments[numberArgument - 1] = fabs(vectorArguments[numberArgument - 1]);

			if (isSameSign)
			{
				minmodVal = (*std::min_element(absVectorArguments.begin(), absVectorArguments.end())) * limiter::mathForLimiter::getSignOfDouble(vectorArguments[0]);
			}
			else
			{
				minmodVal = 0.0;
			}
			return minmodVal;
		}

		double modifiedMinmod(std::vector<double> vectorArguments, double condition)
		{
			double minmodVal(0.0);
			if (fabs(vectorArguments[0])<=condition)
			{
				minmodVal = vectorArguments[0];
			}
			else
			{
				minmodVal = limiter::mathForLimiter::minmod(vectorArguments);
			}
			return minmodVal;
		}

		int getSignOfDouble(double input)
		{
			int sign(0);
			if (input < 0)
			{
				sign = -1;
			}
            else if (input == 0.0)
			{
                sign = 0;
			}
			else
			{
				sign = 1;
			}
			return sign;
		}

		double calM(int element, int valType)
		{
			int neighborElemId(0);
			double M(0.0);
			for (int i = 0; i < auxUlti::checkType(element); i++)
			{
                neighborElemId = meshVar::esuel[element][i];
				if (neighborElemId<0)
				{
					neighborElemId = element;
				}
				switch (valType)
				{
				case 1:
				{
					M += rho[neighborElemId][0] - rho[element][0];
				}
				break;
				case 2:
				{
					M += rhou[neighborElemId][0] - rhou[element][0];
				}
				break;
				case 3:
				{
					M += rhov[neighborElemId][0] - rhov[element][0];
				}
				break;
				case 4:
				{
					M += rhoE[neighborElemId][0] - rhoE[element][0];
				}
				break;
				default:
					break;
				}
			}
			M = fabs(M);
			return M;
		}

		void getNeighborElements()
		{
			for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
			{
				for (int i = 0; i < auxUlti::checkType(nelem); i++)
				{
					meshVar::neighboringElements[nelem][i] = auxUlti::getNeighborElement(nelem, auxUlti::getEdgeHasInputOrderOfElement(nelem, i));
				}
			}
		}
	}
}

namespace troubledCellDetection
{
	bool checkTroubledCell(std::vector<double> InputVector, double condition)
	{
		bool needLimit(true);
        if (fabs(limiter::mathForLimiter::modifiedMinmod(InputVector,condition) - InputVector[0]) > 1e-10)
		{
			needLimit = false;
		}
		else
		{
			needLimit = true;
		}
		return needLimit;
	}

    bool checkTroubleCell_massDiff(int elem)
    {
        //Dung discontinuity dectector (Persson & Paraire)
        bool troubleCell(false);
        std::vector<std::vector<double>>
            SnumGs(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
            SdenGs(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
        double a(0.0), b(0.0), S(0.0);

        for (int na = 0; na <= mathVar::nGauss; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss; nb++)
            {
                std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                math::basisFc(a, b, auxUlti::checkType(elem));  //mathVar::B is changed after this command line is excuted

                SnumGs[na][nb] = pow(mathVar::B[mathVar::orderElem]*rho[elem][mathVar::orderElem],2);
                SdenGs[na][nb]=0.0;
                for (int iorder=0; iorder<mathVar::orderElem; iorder++)
                {
                    SdenGs[na][nb]+=pow(mathVar::B[iorder]*rho[elem][iorder],2);
                }
                SdenGs[na][nb]+=SnumGs[na][nb];
            }
        }

        S=math::volumeInte(SnumGs,elem)/math::volumeInte(SdenGs,elem);
        //Check condition
        if (S>1/pow(mathVar::orderElem+1,4))
        {
            troubleCell=true;
        }
        return troubleCell;
    }
}
