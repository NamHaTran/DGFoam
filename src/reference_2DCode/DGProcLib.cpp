#include "DGProcLib.h"
#include "DGMath.h"
#include "VarDeclaration.h"
#include "DGPostProcessLib.h"
#include "DGAuxUltilitiesLib.h"
#include "dynamicVarDeclaration.h"
#include <tuple>
#include "./boundaryConditions/DGBCsLib.h"
#include <algorithm>
#include "DGIOLib.h"
#include <iostream>
#include <math.h>
#include <mpi.h>
#include "./boundaryConditions/BCSupportFncs.h"
#include "./boundaryConditions/DGBCsLib.h"

#include "./parallelFunctions/GaussPointData.h"
#include "./parallelFunctions/parallelVariables.h"

//Debug
#include "debuggingFuncs.h"

//Limiter
#include "./limiters/limiterController.h"
#include "./limiters/mathFunctions.h"

//Extended Navier-Stokes-Fourier model
#include "./extNSFEqns/FranzDurst/DurstModel.h"
#include "./extNSFEqns/FranzDurst/boundaryConditions.h"

//Non equilibrium BCs
#include "./boundaryConditions/customBCs/nonEquilibriumBCs/nonEqmBCs_GenFuncs.h"

//Time varying BCs
#include "./boundaryConditions/BCSupportFncs.h"
#include "./boundaryConditions/customBCs/timeVaryingBCs/timeVaryingBCs_GenFuncs.h"

namespace meshParam
{
    /**
     * @brief Function generates base Gauss points vectors, which is used to generate Gauss parameters Array for volume Gauss and surface Gauss.
     */
    void genBasedGaussPtsVectors()
    {
        math::getGaussCoeffs(mathVar::xGaussVol,mathVar::wGaussVol,mathVar::nGauss2D);
        math::getGaussLobattoCoeffs(mathVar::xGaussLobattoVol,mathVar::wGaussLobattoVol,mathVar::nGauss2D);  //run GaussLobatto for applying limiter

        math::getGaussCoeffs(mathVar::xGaussSur,mathVar::wGaussSur,mathVar::nGauss1D);
        math::getGaussLobattoCoeffs(mathVar::xGaussLobattoSur,mathVar::wGaussLobattoSur,mathVar::nGauss1D); //run GaussLobatto for applying limiter
    }

    /**
     * @brief Function generates Gauss parameters arrays of all volume Gauss points
     */
	void GaussParam()
	{
        for (int na = 0; na <= mathVar::nGauss2D; na++)  //nGauss is started from 0
		{
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
			{
                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                mathVar::GaussPts[nanb][0] = mathVar::xGaussVol[na];
                mathVar::GaussPts[nanb][1] = mathVar::xGaussVol[nb];
                mathVar::GaussLobattoPts[nanb][0] = mathVar::xGaussLobattoVol[na];
                mathVar::GaussLobattoPts[nanb][1] = mathVar::xGaussLobattoVol[nb];

                mathVar::wGaussPts[nanb][0] = mathVar::wGaussVol[na];
                mathVar::wGaussPts[nanb][1] = mathVar::wGaussVol[nb];
                mathVar::wGaussLobattoPts[nanb][0] = mathVar::wGaussLobattoVol[na];
                mathVar::wGaussLobattoPts[nanb][1] = mathVar::wGaussLobattoVol[nb];
			}
		}
	}

    /**
     * @brief Function calculates Jacobi of all volume Gauss points
     */
	void JacobianParam()
	{
		/*2D Jacobi*/
		double a(0.0), b(0.0);
        for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
            for (int na = 0; na <= mathVar::nGauss2D; na++)  //nGauss is started from 0
			{
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
				{
                    //a = mathVar::GaussPts[na][nb][0];
                    //b = mathVar::GaussPts[na][nb][1];
                    std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                    //meshVar::J2D[ielem][na][nb] = math::J2DCal(ielem, a, b);
                    meshVar::J2D[ielem][calcArrId(na,nb,mathVar::nGauss2D+1)] = math::J2DCal(ielem, a, b);
				}
			}
		}
	}

    /**
     * @brief Function generates array of basis function value at all volume Gauss points at all order (shape mode)
     */
	void basisFcParam()
	{
		double a(0.0), b(0.0);
        for (int na = 0; na <= mathVar::nGauss2D; na++)
		{
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
            {
                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

                //Triangle
                std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                math::basisFc(a, b, 3);
                math::dBasisFc(a, b, 3);
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    mathVar::BPts_Tri[order][nanb] = mathVar::B[order];

                    mathVar::dBaPts_Tri[order][nanb] = mathVar::dBa[order];
                    mathVar::dBbPts_Tri[order][nanb] = mathVar::dBb[order];
                }

                //Quadrilateral
                math::basisFc(a, b, 4);
                math::dBasisFc(a, b, 4);
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    mathVar::BPts_Quad[order][nanb] = mathVar::B[order];

                    mathVar::dBaPts_Quad[order][nanb] = mathVar::dBa[order];
                    mathVar::dBbPts_Quad[order][nanb] = mathVar::dBb[order];
                }
			}
		}
	}

	void derivCoordinates()
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double a(0.0), b(0.0);
        for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
            for (int na = 0; na <= mathVar::nGauss2D; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
				{
                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

                    /*
					a = mathVar::GaussPts[na][nb][0];
					b = mathVar::GaussPts[na][nb][1];
					std::tie(dxa, dxb, dya, dyb) = math::Calc_dxydab(ielem, a, b);
					meshVar::dxa[ielem][na][nb] = dxa;
					meshVar::dxb[ielem][na][nb] = dxb;
					meshVar::dya[ielem][na][nb] = dya;
					meshVar::dyb[ielem][na][nb] = dyb;
                    */
                    std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                    std::tie(dxa, dxb, dya, dyb) = math::Calc_dxydab(ielem, a, b);
                    meshVar::dxa[ielem][nanb] = dxa;
                    meshVar::dxb[ielem][nanb] = dxb;
                    meshVar::dya[ielem][nanb] = dya;
                    meshVar::dyb[ielem][nanb] = dyb;
				}
			}
		}
	}

    /**
     * @brief Function calculates cell metrics:
     *
     * - Geometric Center
     * - Cell size. For Tri cell: Radius of Inscribed Circle of element. For Quad cell: longest diagonal line.
     * - Local Cell size. This size is used for pAdaptive limiter.
     */
	void calcCellMetrics()
	{
		int elemType(0);
		double xCG(0.0), yCG(0.0), cellArea(0.0);
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			elemType = auxUlti::checkType(nelem);
			std::vector<double> xCoor(elemType, 0.0),
				yCoor(elemType, 0.0),
				xSubTriCoor(elemType, 0.0),
				ySubTriCoor(elemType, 0.0);
            for (int i = 0; i < elemType; i++)
			{
				std::tie(xCoor[i], yCoor[i]) = auxUlti::getElemCornerCoord(nelem, i);
			}
            //1. Calculate cell size
            if (elemType==3)
            {
                std::tie(meshVar::cellSize[nelem], cellArea) = math::geometricOp::calROfInscribedCircleOfTriElement(xCoor, yCoor);
            }
            else if (elemType==4)
            {
                std::tie(meshVar::cellSize[nelem], cellArea) = math::geometricOp::calSizeOfQuadElement(xCoor, yCoor);
            }

			//2. Compute geometric center of cell
			std::tie(xCG, yCG) = math::geometricOp::calcGeoCenter(xCoor, yCoor, elemType);

			//3. Compute centroid of cell
			if (elemType == 3)
			{
				meshVar::geoCenter[nelem][0] = xCG;
				meshVar::geoCenter[nelem][1] = yCG;
			}
			else  //quad element
			{
				//3. Compute geometric center of sub-triangles which creates polygon
				std::tie(meshVar::geoCenter[nelem][0], meshVar::geoCenter[nelem][1]) = math::geometricOp::calcQuadCentroid(nelem, xCG, yCG, cellArea);
			}

			//4. Calculate local cell size
			meshVar::localCellSize[nelem] = math::geometricOp::calLocalCellSize(nelem, cellArea);
		}
	}

    /**
     * @brief Function calculates length of edge.
     *
     * This length is also 1D Jacobian of edge.
     */
    void calcEdgeLength()
    {
        //Jacobi 1D = edgeLength/2
        int pt1(-1), pt2(-1);
        for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
        {
            pt1=meshVar::inpoed[iedge][0];
            pt2=meshVar::inpoed[iedge][1];
            meshVar::J1D[iedge]=0.5*math::geometricOp::calDistBetween2Points(meshVar::Points[pt1][0], meshVar::Points[pt1][1], meshVar::Points[pt2][0], meshVar::Points[pt2][1]);
        }
    }

    /**
     * @brief Function calculates coefficients of stiff matrix
     */
	void calcStiffMatrixCoeffs()
	{
        std::vector<std::vector<double>>StiffMatrix(mathVar::orderElem + 1,std::vector<double>(mathVar::orderElem + 1, 0.0));
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				stiffMatrixCoeffs[nelem][iorder] = process::calculateStiffMatrixElement(nelem, iorder, iorder);
                //for (int iorder2 = 0; iorder2 <= mathVar::orderElem; ++iorder2) {
                    //StiffMatrix[iorder][iorder2] = process::calculateStiffMatrixElement(nelem, iorder, iorder2);
                //}
			}
		}
	}

    /**
     * @brief Function calculates distance from cell center to boundary edge
     */
    void calcDistanceFromCenterToBCEdge()
    {
        int globleEdge(0), element(0), tempE;
        for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
        {
            globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
            //int bcType(auxUlti::checkBCTypeOfEdge(globleEdge));
            //if (bcType==meshVar::BCTypeID::wall) //type wall
            //{
                std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(globleEdge);
                meshVar::distanceFromCentroidToBCEdge[ilocalEdge]=math::geometricOp::calcDistanceFromCenterToEdge(element,globleEdge);
            //}
        }
    }

    /**
     * @brief Function finds normal projection of cell center to boundary edge
     */
    void findNormProjectionOfCenterToBCEdge()
    {
        int globleEdge(0), element(0);
        double xC, yC, x1, y1, x2, y2, xN, yN, aN, bN;
        for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
        {
            globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
            int bcType(auxUlti::checkBCTypeOfEdge(globleEdge));
            if (bcType==meshVar::BCTypeID::wall) //type wall
            {
                std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(globleEdge);

                xC=meshVar::geoCenter[element][0];
                yC=meshVar::geoCenter[element][1];
                //Coordinates of 2 vertices of edge
                x1=meshVar::Points[meshVar::inpoed[globleEdge][0]][0];
                y1=meshVar::Points[meshVar::inpoed[globleEdge][0]][1];
                x2=meshVar::Points[meshVar::inpoed[globleEdge][1]][0];
                y2=meshVar::Points[meshVar::inpoed[globleEdge][1]][1];

                std::tie(xN,yN)=math::geometricOp::findNormalProjectionOfPointToEdge(xC,yC,x1,y1,x2,y2);
                std::tie(aN, bN) = math::inverseMapping(element, xN, yN);

                meshVar::normProjectionOfCenterToBCEdge_realSysCoor[ilocalEdge][0]=xN;
                meshVar::normProjectionOfCenterToBCEdge_realSysCoor[ilocalEdge][1]=yN;
                meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[ilocalEdge][0]=aN;
                meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[ilocalEdge][1]=bN;
            }
        }
    }
}

namespace process
{
    /**
     * @brief Function prepares case for parallel running by doing below steps:
     *
     * - Load saved results if loadSavedCase = true
     * - Setup initial conditions if loadSavedCase = false
     * - Apply limiter
     * - Send/Receive data between cores
     */
    void prepareCaseForRunning()
    {
        if (systemVar::loadSavedCase)
        {
            meshParam::calcStiffMatrixCoeffs();
            if (systemVar::currentProc==0){
            std::cout << "Loading case...\n" << std::endl;}

            IO::loadCase(systemVar::readWriteMode);

            if (
                    /* Neu khong co cac file TSurface, uSurface, vSurface thi phai
                     * set initial value cua cac gia tri trong SurfaceBCFields
                    */
                    !controlFlag::fileAvailFlags::fileTSurface
                    ||!controlFlag::fileAvailFlags::fileuSurface
                    ||!controlFlag::fileAvailFlags::filevSurface
                    /* Hoac neu khong co dieu kien time-varying cung phai set initial
                     * value
                    */
                    ||!auxUlti::checkTimeVaryingBCAvailable())
            {
                process::setIniSurfaceBCValues();
            }
        }
        else
        {
            /*SET INITIAL VALUES*/
            if (systemVar::initializedOrNot == false)
            {
                meshParam::calcStiffMatrixCoeffs();
                process::setIniVolumeValues();
                process::setIniSurfaceBCValues();
            }
        }

        BCSupportFncs::calcBCValuesAtUnfixedValueBCs();

        if (systemVar::currentProc==0)
        {
            std::cout << " \n" << "Simulation is started\n";
        }

        //SETUP FOR LIMITER
        /* Khong apply limiter o buoc nay*/
        limiter::mathForLimiter::getNeighborElements();
        limitVal::numOfLimitCell = 0;

        // limiter::limiter_1OutterStep(); -> cause crash becase of very large u and v values
        // Nguyen nhan vi ham mass diffusion limiter dung chung BR1Vars::rhoX va BR1Vars::rhoY voi ham giai pt phu S = mu(div(U))
        // Nen xay ra conflict la mass diffusion limiter luu div(rho) va khi giai pt phu thi luu mu*div(rho)

        limiter::limiter_1InnerStep();
    }

    /*Ham lam cac thu tuc truoc khi ket thuc 1 iteration*/
    void finalizeIteration()
    {
        if (systemVar::savingCout >= systemVar::wrtI)
        {
            std::cout << "Saving case...\n" << std::endl;
            IO::saveCase();
            std::cout << "Exporting data to Tecplot...\n" << std::endl;
            DG2Tecplot::exportCellCenteredData(systemVar::iterCount);
            systemVar::savingCout = 0;
        }

        if (systemVar::loadConstCount == 10)
        {
            IO::loadSettingFiles::loadConstantsWhileRunning();
            systemVar::loadConstCount = 0;
        }
        systemVar::loadConstCount++;

        //set lai cac bien flag
        systemVar::firstIter=false;
        if (mathVar::solveTFailed)
        {
            std::cout<<"Failed to solve T at somewhere!\n";
            mathVar::solveTFailed=false;
        }
    }

    void setIniVolumeValues()
	{
		iniValues::rhoIni = iniValues::pIni / (material::R*iniValues::TIni);
		iniValues::eIni = material::Cv*iniValues::TIni; //+0.5*(pow(iniValues::uIni, 2) + pow(iniValues::vIni, 2) + pow(iniValues::wIni, 2));
		iniValues::muIni = math::CalcVisCoef(iniValues::TIni);
		std::vector<double> iniRho(mathVar::orderElem + 1),
			iniRhou(mathVar::orderElem + 1),
			iniRhov(mathVar::orderElem + 1),
			iniRhoE(mathVar::orderElem + 1);

        for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
		{
			iniRho = process::calcIniValues(iniValues::rhoIni, nelement);
			iniRhou = process::calcIniValues(iniValues::rhoIni*iniValues::uIni, nelement);
			iniRhov = process::calcIniValues(iniValues::rhoIni*iniValues::vIni, nelement);
			iniRhoE = process::calcIniValues(iniValues::rhoIni*(iniValues::eIni + 0.5*(pow(iniValues::uIni, 2) + pow(iniValues::vIni, 2))), nelement);
            for (int i = 0; i <= mathVar::orderElem; i++)
			{
				rho[nelement][i] = iniRho[i];
				rhou[nelement][i] = iniRhou[i];
				rhov[nelement][i] = iniRhov[i];
				rhoE[nelement][i] = iniRhoE[i];
            }
        }

        //Calculate limit of rhoE
        //limitVal::rhoEUp = material::Cv*limitVal::TUp*limitVal::rhoUp;
        //limitVal::rhoEDwn = material::Cv*limitVal::TDwn*limitVal::rhoDwn;
	}

    void setIniSurfaceBCValues()
    {
        int BCGrp(0), globleEdgeId(0);
        for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
        {
            globleEdgeId=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
            BCGrp=auxUlti::getGrpOfEdge(globleEdgeId);

            for (int nG=0; nG<=mathVar::nGauss1D; nG++)
            {
                SurfaceBCFields::TBc[ilocalEdge][nG]=bcValues::TBCFixed[BCGrp-1];
                SurfaceBCFields::uBc[ilocalEdge][nG]=bcValues::uBCFixed[BCGrp-1];
                SurfaceBCFields::vBc[ilocalEdge][nG]=bcValues::vBCFixed[BCGrp-1];
                SurfaceBCFields::pBc[ilocalEdge][nG]=bcValues::pBCFixed[BCGrp-1];
            }
        }
    }

	std::vector<double> calcIniValues(double iniVal, int element)
	{
		std::vector<double> RHS(mathVar::orderElem + 1, 0.0);
		std::vector<double> iniVector(mathVar::orderElem + 1, 0.0);

		RHS = process::calcIniValuesRHS(element, iniVal);
        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
		{
			iniVector[iorder] = RHS[iorder] / stiffMatrixCoeffs[element][iorder];
		}
		return iniVector;
	}

	std::vector<double> calcIniValuesRHS(int element, double iniVal)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
        std::vector<std::vector<double>> matrix(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));
		double a(0.0), b(0.0);
        for (int order = 0; order <= mathVar::orderElem; order++)
		{
            for (int na = 0; na <= mathVar::nGauss2D; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					math::basisFc(a, b, auxUlti::checkType(element));
					matrix[na][nb] = mathVar::B[order];
				}
			}
			Out[order] = math::volumeInte(matrix, element)*iniVal;
		}
		return Out;
    }

    /**
     * @brief Function calculates conservative variables at all volume Gauss points
     */
	void calcVolumeGaussValues()
	{
        /* Ham tinh gia tri U tai cac diem Gauss ben trong cell
         * Ket qua duoc luu vao volumeFields
        */
		double a(0.0), b(0.0);
        std::vector<double>UVol(4,0.0);

        if (mathVar::orderElem>0)
        {
            for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
            {
                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

                        std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                        UVol=math::pointUVars(nelem,a,b);

                        volumeFields::rhoVolGauss[nelem][nanb] = UVol[0];
                        volumeFields::rhouVolGauss[nelem][nanb] = UVol[1];
                        volumeFields::rhovVolGauss[nelem][nanb] = UVol[2];
                        volumeFields::rhoEVolGauss[nelem][nanb] = UVol[3];
                    }
                }
            }
        }
        /* Neu order = 0 thi khong can goi ham pointUVars*/
        else
        {
            for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
            {
                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

                        volumeFields::rhoVolGauss[nelem][nanb] = rho[nelem][0];
                        volumeFields::rhouVolGauss[nelem][nanb] = rhou[nelem][0];
                        volumeFields::rhovVolGauss[nelem][nanb] = rhov[nelem][0];
                        volumeFields::rhoEVolGauss[nelem][nanb] = rhoE[nelem][0];
                    }
                }
            }
        }
	}

    /**
     * @brief Function calculates Rho at all volume Gauss points.
     *
     * Values are saved in volumeFields::rhoVolGauss.
     */
    void calcVolumeGaussRho()
    {
        /* Ham tinh gia tri U tai cac diem Gauss ben trong cell
         * Ket qua duoc luu vao volumeFields
        */
        double a(0.0), b(0.0);

        if (mathVar::orderElem>0)
        {
            for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
            {
                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

                        std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                        volumeFields::rhoVolGauss[nelem][nanb] = math::pointValue(nelem,a,b,1,2);
                    }
                }
            }
        }
        /* Neu order = 0 thi khong can goi ham pointUVars*/
        else
        {
            for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
            {
                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                        volumeFields::rhoVolGauss[nelem][nanb] = rho[nelem][0];
                    }
                }
            }
        }
    }

    void calcTGauss()
    {
        //Volume
        double a(0.0), b(0.0);
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
        {
            for (int na = 0; na <= mathVar::nGauss2D; na++)
            {
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                {
                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

                    std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                    volumeFields::T[nelem][nanb] = math::pointValue(nelem, a, b, 6, 1);
                }
            }
        }

        //Surface
        int masterCell(-1), slaveCell(-1), bcGrp(0);
        double rhoMaster(0.0), rhouMaster(0.0), rhovMaster(0.0), rhoEMaster(0.0),
                rhoSlave(0.0), rhouSlave(0.0), rhovSlave(0.0), rhoESlave(0.0);
        for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
        {
            bcGrp = auxUlti::getBCType(iedge);
            if (bcGrp == 0)
            {
                std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    std::tie(rhoMaster,rhoSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,1);
                    std::tie(rhouMaster,rhouSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,2);
                    std::tie(rhovMaster,rhovSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,3);
                    std::tie(rhoEMaster,rhoESlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,4);
                    surfaceFields::T[iedge][nG]=math::CalcTFromConsvVar(rhoMaster,rhouMaster,rhovMaster,rhoEMaster);
                    surfaceFields::T[iedge][nG + mathVar::nGauss1D + 1]=math::CalcTFromConsvVar(rhoSlave,rhouSlave,rhovSlave,rhoESlave);
                }
            }
            else
            {
                std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    rhoMaster = surfaceFields::rho[iedge][nG];
                    rhouMaster = surfaceFields::rhou[iedge][nG];
                    rhovMaster = surfaceFields::rhov[iedge][nG];
                    rhoEMaster = surfaceFields::rhoE[iedge][nG];

                    surfaceFields::T[iedge][nG]=math::CalcTFromConsvVar(rhoMaster,rhouMaster,rhovMaster,rhoEMaster);
                }
            }
        }

        if (systemVar::parallelMode)
        {
            parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::T,surfaceFields::T);
        }
    }

    void calcValuesAtInterface()
    {
        int masterCell(-1), slaveCell(-1), bcGrp(0);
        double tempMaster(0.0), tempSlave(0.0);

        if (mathVar::orderElem>0)
        {
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
            {
                bcGrp = auxUlti::getBCType(iedge);
                if (bcGrp == 0)
                {
                    std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        //rho
                        std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 1, 2);
                        surfaceFields::rho[iedge][nG] = tempMaster;
                        surfaceFields::rho[iedge][nG + mathVar::nGauss1D + 1] = tempSlave;

                        //rhou
                        std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 2, 2);
                        surfaceFields::rhou[iedge][nG] = tempMaster;
                        surfaceFields::rhou[iedge][nG + mathVar::nGauss1D + 1] = tempSlave;

                        //rhov
                        std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 3, 2);
                        surfaceFields::rhov[iedge][nG] = tempMaster;
                        surfaceFields::rhov[iedge][nG + mathVar::nGauss1D + 1] = tempSlave;

                        //rhoE
                        std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 4, 2);
                        surfaceFields::rhoE[iedge][nG] = tempMaster;
                        surfaceFields::rhoE[iedge][nG + mathVar::nGauss1D + 1] = tempSlave;
                    }
                }
                /*
                else if (bcGrp == 4)
                {
                    std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        //rho
                        surfaceFields::rho[iedge][nG] = rho[masterCell][0];

                        //rhou
                        surfaceFields::rhou[iedge][nG] = rhou[masterCell][0];

                        //rhov
                        surfaceFields::rhov[iedge][nG] = rhov[masterCell][0];

                        //rhoE
                        surfaceFields::rhoE[iedge][nG] = rhoE[masterCell][0];
                    }
                }*/
                else
                {
                    std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        //rho
                        surfaceFields::rho[iedge][nG] = math::plusSideSurfaceValue(iedge, masterCell, nG, 1, 2);

                        //rhou
                        surfaceFields::rhou[iedge][nG] = math::plusSideSurfaceValue(iedge, masterCell, nG, 2, 2);

                        //rhov
                        surfaceFields::rhov[iedge][nG] = math::plusSideSurfaceValue(iedge, masterCell, nG, 3, 2);

                        //rhoE
                        surfaceFields::rhoE[iedge][nG] = math::plusSideSurfaceValue(iedge, masterCell, nG, 4, 2);
                    }
                }
            }
        }
        /*Khong can goi ham internalSurfaceValue khi order = 0*/
        else
        {
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
            {
                bcGrp = auxUlti::getBCType(iedge);
                if (bcGrp == 0)
                {
                    std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        //rho
                        surfaceFields::rho[iedge][nG] = rho[masterCell][0];
                        surfaceFields::rho[iedge][nG + mathVar::nGauss1D + 1] = rho[slaveCell][0];

                        //rhou
                        surfaceFields::rhou[iedge][nG] = rhou[masterCell][0];
                        surfaceFields::rhou[iedge][nG + mathVar::nGauss1D + 1] = rhou[slaveCell][0];

                        //rhov
                        surfaceFields::rhov[iedge][nG] = rhov[masterCell][0];
                        surfaceFields::rhov[iedge][nG + mathVar::nGauss1D + 1] = rhov[slaveCell][0];

                        //rhoE
                        surfaceFields::rhoE[iedge][nG] = rhoE[masterCell][0];
                        surfaceFields::rhoE[iedge][nG + mathVar::nGauss1D + 1] = rhoE[slaveCell][0];
                    }
                }
                else
                {
                    std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        //rho
                        surfaceFields::rho[iedge][nG] = rho[masterCell][0];

                        //rhou
                        surfaceFields::rhou[iedge][nG] = rhou[masterCell][0];

                        //rhov
                        surfaceFields::rhov[iedge][nG] = rhov[masterCell][0];

                        //rhoE
                        surfaceFields::rhoE[iedge][nG] = rhoE[masterCell][0];
                    }
                }
            }
        }

        if (systemVar::parallelMode)
        {
            parallelFuncs_GaussPt::synchSurfGaussU();
        }
    }

    /**
     * @brief Function calculates Rho at all surface Gauss points.
     *
     * Values are saved in surfaceFields::rho.\n
     * - Master cell side: [...][nG]\n
     * - Slave cell side: [...][nG + mathVar::nGauss + 1]
     */
    void calcRhoAtInterface()
    {
        int masterCell(-1), slaveCell(-1), bcGrp(0);

        if (mathVar::orderElem>0)
        {
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
            {
                bcGrp = auxUlti::getBCType(iedge);
                if (bcGrp == 0)
                {
                    std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        std::tie(surfaceFields::rho[iedge][nG], surfaceFields::rho[iedge][nG + mathVar::nGauss1D + 1])
                                = math::internalSurfaceValue(iedge, masterCell, nG, 1, 2);
                    }
                }
                /*
                else if (bcGrp == 4)
                {
                    std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        surfaceFields::rho[iedge][nG] = math::plusSideSurfaceValue(iedge, masterCell, nG, 1, 2);
                    }
                }*/
                else
                {
                    std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        surfaceFields::rho[iedge][nG] = math::plusSideSurfaceValue(iedge, masterCell, nG, 1, 2);
                    }
                }
            }
        }
        /*Khong can goi ham internalSurfaceValue khi order = 0*/
        else
        {
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
            {
                bcGrp = auxUlti::getBCType(iedge);
                if (bcGrp == 0)
                {
                    std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        //rho
                        surfaceFields::rho[iedge][nG] = rho[masterCell][0];
                        surfaceFields::rho[iedge][nG + mathVar::nGauss1D + 1] = rho[slaveCell][0];
                    }
                }
                else
                {
                    std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        //rho
                        surfaceFields::rho[iedge][nG] = rho[masterCell][0];
                    }
                }
            }
        }

        if (systemVar::parallelMode)
        {
            parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::rho,surfaceFields::rho);
        }
    }

	namespace auxEq
    {
    void updateSurfaceFieldsAtBC()
    {

        //Ham su dung khi flow la inviscid, co chuc nang update surface fields tai cac BC edge de su dung cho phan giai phuong trinh chinh
        for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
        {
            int elemType(auxUlti::checkType(nelement)), edgeName(0);
            int faceBcType(0);

            for (int nface = 0; nface < elemType; nface++)
            {
                edgeName = meshVar::inedel[nelement][nface];
                faceBcType = auxUlti::getBCType(edgeName);
                if (faceBcType != 0)  //boundary edge
                {
                    for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
                    {
                        //Chi chay ham auxEqBCsImplement de tinh toan va cap nhat U trong surfaceFields tai BC
                        auxEqBCsImplement(nelement, edgeName, nGauss);
                    }
                }
            }
        }
    }

    void solveAuxEquation()
    {
        if (systemVar::auxVariables==1)
        {
            //std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
            std::vector<double> rhoRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhoRHSTermOyDir(mathVar::orderElem + 1, 0.0),
                rhouRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhouRHSTermOyDir(mathVar::orderElem + 1, 0.0),
                rhovRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhovRHSTermOyDir(mathVar::orderElem + 1, 0.0),
                rhoERHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhoERHSTermOyDir(mathVar::orderElem + 1, 0.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
            {
                //1) Calculate Right hand side terms
                process::auxEq::BR1::CalcRHSTerm(nelement, rhoRHSTermOxDir, rhoRHSTermOyDir, rhouRHSTermOxDir, rhouRHSTermOyDir, rhovRHSTermOxDir, rhovRHSTermOyDir, rhoERHSTermOxDir, rhoERHSTermOyDir);

                //2) Solve for auxilary variables, NOTE!!: bien phu la mu*div(U)
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    //Ox direction
                    BR1Vars::rhoX[nelement][iorder] = rhoRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhouX[nelement][iorder] = rhouRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhovX[nelement][iorder] = rhovRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhoEX[nelement][iorder] = rhoERHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];

                    //Oy direction
                    BR1Vars::rhoY[nelement][iorder] = rhoRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhouY[nelement][iorder] = rhouRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhovY[nelement][iorder] = rhovRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhoEY[nelement][iorder] = rhoERHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                }
            }
        }
        else if (systemVar::auxVariables==2)
        {
            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
            {
                //1) Calculate volume divU
                process::auxEq::BR2::calcVolDivU(nelement);

                //2) Calculate suface divU
                process::auxEq::BR2::calcSurDivU(nelement);
            }
        }
    }
    namespace BR2 {
    void calcVolDivU(int element)
    {
        std::vector<std::vector<double>> gaussMX(mathVar::nGauss2D+1,std::vector<double>(mathVar::nGauss2D+1,0.0)),
                gaussMY(mathVar::nGauss2D+1,std::vector<double>(mathVar::nGauss2D+1,0.0));
        std::vector<double>rhoC(mathVar::orderElem+1,0.0),
                rhouC(mathVar::orderElem+1,0.0),
                rhovC(mathVar::orderElem+1,0.0),
                rhoEC(mathVar::orderElem+1,0.0);
        int elemType(auxUlti::checkType(element));
        double Bn(0.0),dBmX(0.0),dBmY(0.0),
                GammaRhoX(0.0), GammaRhoY(0.0),
                GammaRhouX(0.0), GammaRhouY(0.0),
                GammaRhovX(0.0), GammaRhovY(0.0),
                GammaRhoEX(0.0), GammaRhoEY(0.0),
                GammaX(0.0), GammaY(0.0), theta1(1.0), theta2(1.0);

        rhoC=auxUlti::getElementConserValuesOfOrder(element,1);
        rhouC=auxUlti::getElementConserValuesOfOrder(element,2);
        rhovC=auxUlti::getElementConserValuesOfOrder(element,3);
        rhoEC=auxUlti::getElementConserValuesOfOrder(element,4);
        for (int n=0;n<=mathVar::orderElem;n++) {
            GammaRhoX=0; GammaRhoY=0;
            GammaRhouX=0; GammaRhouY=0;
            GammaRhovX=0; GammaRhovY=0;
            GammaRhoEX=0; GammaRhoEY=0;

            for (int m=0;m<=mathVar::orderElem;m++) {
                if (m==0)
                {
                    theta1=1.0;
                    theta2=1.0;
                }
                else {
                    theta1=theta1Arr[element];
                    theta2=theta2Arr[element];
                }

                if (elemType==4)
                {
                    for (int na=0;na<=mathVar::nGauss2D;na++) {
                        for (int nb=0;nb<=mathVar::nGauss2D;nb++) {
                            int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                            Bn=mathVar::BPts_Quad[n][nanb];
                            dBmX=math::Calc_dBxdBy(element, m, na, nb, 1);
                            dBmY=math::Calc_dBxdBy(element, m, na, nb, 2);
                            gaussMX[na][nb]=Bn*dBmX;
                            gaussMY[na][nb]=Bn*dBmY;
                        }
                    }
                }
                else if (elemType==3)
                {
                    for (int na=0;na<=mathVar::nGauss2D;na++) {
                        for (int nb=0;nb<=mathVar::nGauss2D;nb++) {
                            int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                            Bn=mathVar::BPts_Tri[n][nanb];
                            dBmX=math::Calc_dBxdBy(element, m, na, nb, 1);
                            dBmY=math::Calc_dBxdBy(element, m, na, nb, 2);
                            gaussMX[na][nb]=Bn*dBmX;
                            gaussMY[na][nb]=Bn*dBmY;
                        }
                    }
                }

                GammaX=math::volumeInte(gaussMX,element);
                GammaY=math::volumeInte(gaussMY,element);

                GammaRhoX+=GammaX*rhoC[m]*theta1*theta2;
                GammaRhoY+=GammaY*rhoC[m]*theta1*theta2;

                GammaRhouX+=GammaX*rhouC[m]*theta2;
                GammaRhouY+=GammaY*rhouC[m]*theta2;

                GammaRhovX+=GammaX*rhovC[m]*theta2;
                GammaRhovY+=GammaY*rhovC[m]*theta2;

                GammaRhoEX+=GammaX*rhoEC[m]*theta2;
                GammaRhoEY+=GammaY*rhoEC[m]*theta2;
            }
            double stiffCoeff(stiffMatrixCoeffs[element][n]);
            BR2Vars::rhoXVol[element][n]=GammaRhoX/stiffCoeff;
            BR2Vars::rhoYVol[element][n]=GammaRhoY/stiffCoeff;

            BR2Vars::rhouXVol[element][n]=GammaRhouX/stiffCoeff;
            BR2Vars::rhouYVol[element][n]=GammaRhouY/stiffCoeff;

            BR2Vars::rhovXVol[element][n]=GammaRhovX/stiffCoeff;
            BR2Vars::rhovYVol[element][n]=GammaRhovY/stiffCoeff;

            BR2Vars::rhoEXVol[element][n]=GammaRhoEX/stiffCoeff;
            BR2Vars::rhoEYVol[element][n]=GammaRhoEY/stiffCoeff;
        }
    }

    void calcSurDivU(int element)
    {
        int elemType(auxUlti::checkType(element)), edgeId(0), BCType(-1);
        double a(0.0), b(0.0),BVar(0.0), nx(0.0), ny(0.0), UP(0.0), UM(0.0), BRConst(4.0),
                dRhoXVar(0.0), dRhoYVar(0.0),
                dRhouXVar(0.0), dRhouYVar(0.0),
                dRhovXVar(0.0), dRhovYVar(0.0),
                dRhoEXVar(0.0), dRhoEYVar(0.0);
        bool isMaster;
        std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0));
        std::vector<double> vectorRhoX(mathVar::nGauss1D+1,0.0),
                vectorRhoY(mathVar::nGauss1D+1,0.0),
                vectorRhouX(mathVar::nGauss1D+1,0.0),
                vectorRhouY(mathVar::nGauss1D+1,0.0),
                vectorRhovX(mathVar::nGauss1D+1,0.0),
                vectorRhovY(mathVar::nGauss1D+1,0.0),
                vectorRhoEX(mathVar::nGauss1D+1,0.0),
                vectorRhoEY(mathVar::nGauss1D+1,0.0),

                vectordRhoXVol_org(mathVar::orderElem+1,0.0),
                vectordRhoYVol_org(mathVar::orderElem+1,0.0),
                vectordRhouXVol_org(mathVar::orderElem+1,0.0),
                vectordRhouYVol_org(mathVar::orderElem+1,0.0),
                vectordRhovXVol_org(mathVar::orderElem+1,0.0),
                vectordRhovYVol_org(mathVar::orderElem+1,0.0),
                vectordRhoEXVol_org(mathVar::orderElem+1,0.0),
                vectordRhoEYVol_org(mathVar::orderElem+1,0.0);

        for (int iorder=0;iorder<=mathVar::orderElem;iorder++) {
            vectordRhoXVol_org[iorder]=BR2Vars::rhoXVol[element][iorder];
            vectordRhoYVol_org[iorder]=BR2Vars::rhoYVol[element][iorder];
            vectordRhouXVol_org[iorder]=BR2Vars::rhouXVol[element][iorder];
            vectordRhouYVol_org[iorder]=BR2Vars::rhouYVol[element][iorder];
            vectordRhovXVol_org[iorder]=BR2Vars::rhovXVol[element][iorder];
            vectordRhovYVol_org[iorder]=BR2Vars::rhovYVol[element][iorder];
            vectordRhoEXVol_org[iorder]=BR2Vars::rhoEXVol[element][iorder];
            vectordRhoEYVol_org[iorder]=BR2Vars::rhoEYVol[element][iorder];
        }

        for (int iedge=0;iedge<elemType;iedge++) {
            edgeId=meshVar::inedel[element][iedge];
            BCType=auxUlti::getBCType(edgeId);
            isMaster=auxUlti::checkMaster(element,edgeId);
            nx=auxUlti::getNormVectorComp(element,edgeId,1);
            ny=auxUlti::getNormVectorComp(element,edgeId,2);

            if (BCType!=0)
            {
                for (int nG=0;nG<=mathVar::nGauss1D;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                    gaussVector = auxEqBCsImplement(element, edgeId, nG);

                    //Rho
                    vectorRhoX[nG]=0.5*(gaussVector[0][1]-gaussVector[0][0])*nx;
                    vectorRhoY[nG]=0.5*(gaussVector[0][1]-gaussVector[0][0])*ny;
                }
            }
            else {
                for (int nG=0;nG<=mathVar::nGauss1D;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);

                    //Rho
                    std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,1);
                    vectorRhoX[nG]=0.5*(UM-UP)*nx;
                    vectorRhoY[nG]=0.5*(UM-UP)*ny;
                }
            }

            if (BCType!=0)
            {
                for (int nG=0;nG<=mathVar::nGauss1D;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                    gaussVector = auxEqBCsImplement(element, edgeId, nG);

                    //Rhou
                    vectorRhouX[nG]=0.5*(gaussVector[1][1]-gaussVector[1][0])*nx;
                    vectorRhouY[nG]=0.5*(gaussVector[1][1]-gaussVector[1][0])*ny;

                    //Rhov
                    vectorRhovX[nG]=0.5*(gaussVector[2][1]-gaussVector[2][0])*nx;
                    vectorRhovY[nG]=0.5*(gaussVector[2][1]-gaussVector[2][0])*ny;

                    //RhoE
                    vectorRhoEX[nG]=0.5*(gaussVector[3][1]-gaussVector[3][0])*nx;
                    vectorRhoEY[nG]=0.5*(gaussVector[3][1]-gaussVector[3][0])*ny;
                }
            }
            else {
                for (int nG=0;nG<=mathVar::nGauss1D;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);

                    //Rhou
                    std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,2);
                    vectorRhouX[nG]=0.5*(UM-UP)*nx;
                    vectorRhouY[nG]=0.5*(UM-UP)*ny;

                    //Rhov
                    std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,3);
                    vectorRhovX[nG]=0.5*(UM-UP)*nx;
                    vectorRhovY[nG]=0.5*(UM-UP)*ny;

                    //RhoE
                    std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,4);
                    vectorRhoEX[nG]=0.5*(UM-UP)*nx;
                    vectorRhoEY[nG]=0.5*(UM-UP)*ny;
                }
            }

            for (int iorder=0;iorder<=mathVar::orderElem;iorder++) {
                double stiffCoeff(stiffMatrixCoeffs[element][iorder]);
                for (int nG=0;nG<=mathVar::nGauss1D;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                    math::basisFc(a,b,elemType);
                    BVar=mathVar::B[iorder];
                    //Rho
                    vectorRhoX[nG]*=BVar;
                    vectorRhoY[nG]*=BVar;

                    //Rhou
                    vectorRhouX[nG]*=BVar;
                    vectorRhouY[nG]*=BVar;

                    //Rhov
                    vectorRhovX[nG]*=BVar;
                    vectorRhovY[nG]*=BVar;

                    //RhoE
                    vectorRhoEX[nG]*=BVar;
                    vectorRhoEY[nG]*=BVar;
                }

                dRhoXVar=math::surfaceInte(vectorRhoX,edgeId)/stiffCoeff;
                dRhoYVar=math::surfaceInte(vectorRhoY,edgeId)/stiffCoeff;

                dRhouXVar=math::surfaceInte(vectorRhouX,edgeId)/stiffCoeff;
                dRhouYVar=math::surfaceInte(vectorRhouY,edgeId)/stiffCoeff;

                dRhovXVar=math::surfaceInte(vectorRhovX,edgeId)/stiffCoeff;
                dRhovYVar=math::surfaceInte(vectorRhovY,edgeId)/stiffCoeff;

                dRhoEXVar=math::surfaceInte(vectorRhoEX,edgeId)/stiffCoeff;
                dRhoEYVar=math::surfaceInte(vectorRhoEY,edgeId)/stiffCoeff;

                if (isMaster)
                {
                    BR2Vars::rhoXSurMaster[edgeId][iorder]=BRConst*dRhoXVar+vectordRhoXVol_org[iorder];
                    BR2Vars::rhoYSurMaster[edgeId][iorder]=BRConst*dRhoYVar+vectordRhoYVol_org[iorder];

                    BR2Vars::rhouXSurMaster[edgeId][iorder]=BRConst*dRhouXVar+vectordRhouXVol_org[iorder];
                    BR2Vars::rhouYSurMaster[edgeId][iorder]=BRConst*dRhouYVar+vectordRhouYVol_org[iorder];

                    BR2Vars::rhovXSurMaster[edgeId][iorder]=BRConst*dRhovXVar+vectordRhovXVol_org[iorder];
                    BR2Vars::rhovYSurMaster[edgeId][iorder]=BRConst*dRhovYVar+vectordRhovYVol_org[iorder];

                    BR2Vars::rhoEXSurMaster[edgeId][iorder]=BRConst*dRhoEXVar+vectordRhoEXVol_org[iorder];
                    BR2Vars::rhoEYSurMaster[edgeId][iorder]=BRConst*dRhoEYVar+vectordRhoEYVol_org[iorder];
                }
                else {
                    BR2Vars::rhoXSurSlave[edgeId][iorder]=BRConst*dRhoXVar+vectordRhoXVol_org[iorder];
                    BR2Vars::rhoYSurSlave[edgeId][iorder]=BRConst*dRhoYVar+vectordRhoYVol_org[iorder];

                    BR2Vars::rhouXSurSlave[edgeId][iorder]=BRConst*dRhouXVar+vectordRhouXVol_org[iorder];
                    BR2Vars::rhouYSurSlave[edgeId][iorder]=BRConst*dRhouYVar+vectordRhouYVol_org[iorder];

                    BR2Vars::rhovXSurSlave[edgeId][iorder]=BRConst*dRhovXVar+vectordRhovXVol_org[iorder];
                    BR2Vars::rhovYSurSlave[edgeId][iorder]=BRConst*dRhovYVar+vectordRhovYVol_org[iorder];

                    BR2Vars::rhoEXSurSlave[edgeId][iorder]=BRConst*dRhoEXVar+vectordRhoEXVol_org[iorder];
                    BR2Vars::rhoEYSurSlave[edgeId][iorder]=BRConst*dRhoEYVar+vectordRhoEYVol_org[iorder];
                }

                BR2Vars::rhoXVol[element][iorder]+=dRhoXVar;
                BR2Vars::rhoYVol[element][iorder]+=dRhoYVar;

                BR2Vars::rhouXVol[element][iorder]+=dRhouXVar;
                BR2Vars::rhouYVol[element][iorder]+=dRhouYVar;

                BR2Vars::rhovXVol[element][iorder]+=dRhovXVar;
                BR2Vars::rhovYVol[element][iorder]+=dRhovYVar;

                BR2Vars::rhoEXVol[element][iorder]+=dRhoEXVar;
                BR2Vars::rhoEYVol[element][iorder]+=dRhoEYVar;
            }
        }
    }


    }
    namespace BR1 {
    void calcValuesAtInterface()
    {
        //mu included
        int masterCell(-1), slaveCell(-1), bcGrp(0);
        double muMaster(0.0), muSlave(0.0), tempMaster(0.0), tempSlave(0.0),TMaster(0.0), TSlave(0.0);
        for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
        {
            bcGrp = auxUlti::getBCType(iedge);
            if (bcGrp == 0)
            {
                std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    std::tie(TMaster, TSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 6, 1);
                    muMaster=math::CalcVisCoef(TMaster);
                    muSlave=math::CalcVisCoef(TSlave);
                    surfaceFields::T[iedge][nG]=TMaster;
                    surfaceFields::T[iedge][nG+ mathVar::nGauss1D + 1]=TSlave;

                    //rho*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 1);
                    surfaceFields::aux_rho[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rho[iedge][nG + mathVar::nGauss1D + 1] = tempSlave * muSlave;

                    //rhou*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 2);
                    surfaceFields::aux_rhou[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rhou[iedge][nG + mathVar::nGauss1D + 1] = tempSlave * muSlave;

                    //rhov*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 3);
                    surfaceFields::aux_rhov[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rhov[iedge][nG + mathVar::nGauss1D + 1] = tempSlave * muSlave;

                    //rhoE*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 4);
                    surfaceFields::aux_rhoE[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rhoE[iedge][nG + mathVar::nGauss1D + 1] = tempSlave * muSlave;
                }
            }
        }
    }

    std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType)
    {
        //valType = 1(rho*mu), 2(rhou*mu), 3(rhov*mu), 4(rhoE*mu)
        bool isMaster(auxUlti::checkMaster(element, edge));
        int locationPlus(-1), locationMinus(-1);
        double valPlus(0.0), valMinus(0.0);
        if (isMaster)
        {
            locationPlus = nG;
            locationMinus = nG + mathVar::nGauss1D + 1;
        }
        else
        {
            locationPlus = nG + mathVar::nGauss1D + 1;
            locationMinus = nG;
        }

        switch (valType)
        {
        case 1:
        {
            valPlus = surfaceFields::aux_rho[edge][locationPlus];
            valMinus = surfaceFields::aux_rho[edge][locationMinus];
        }
        break;
        case 2:
        {
            valPlus = surfaceFields::aux_rhou[edge][locationPlus];
            valMinus = surfaceFields::aux_rhou[edge][locationMinus];
        }
        break;
        case 3:
        {
            valPlus = surfaceFields::aux_rhov[edge][locationPlus];
            valMinus = surfaceFields::aux_rhov[edge][locationMinus];
        }
        break;
        case 4:
        {
            valPlus = surfaceFields::aux_rhoE[edge][locationPlus];
            valMinus = surfaceFields::aux_rhoE[edge][locationMinus];
        }
        break;
        default:
            break;
        }
        return std::make_tuple(valPlus, valMinus);
    }

    /*Function calculates right hand side terms of all conservative variables at all order in all directions*/
    void CalcRHSTerm(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy, std::vector<double> &rhouRHSOx, std::vector<double> &rhouRHSOy, std::vector<double> &rhovRHSOx, std::vector<double> &rhovRHSOy, std::vector<double> &rhoERHSOx, std::vector<double> &rhoERHSOy)
    {
        std::vector<double>
            //vectors for volume integrals
            rhoVolIntOx(mathVar::orderElem + 1, 0.0),
            rhoVolIntOy(mathVar::orderElem + 1, 0.0),
            rhouVolIntOx(mathVar::orderElem + 1, 0.0),
            rhouVolIntOy(mathVar::orderElem + 1, 0.0),
            rhovVolIntOx(mathVar::orderElem + 1, 0.0),
            rhovVolIntOy(mathVar::orderElem + 1, 0.0),
            rhoEVolIntOx(mathVar::orderElem + 1, 0.0),
            rhoEVolIntOy(mathVar::orderElem + 1, 0.0),

            //vectors for surface integrals
            rhoSurfIntOx(mathVar::orderElem + 1, 0.0),
            rhoSurfIntOy(mathVar::orderElem + 1, 0.0),
            rhouSurfIntOx(mathVar::orderElem + 1, 0.0),
            rhouSurfIntOy(mathVar::orderElem + 1, 0.0),
            rhovSurfIntOx(mathVar::orderElem + 1, 0.0),
            rhovSurfIntOy(mathVar::orderElem + 1, 0.0),
            rhoESurfIntOx(mathVar::orderElem + 1, 0.0),
            rhoESurfIntOy(mathVar::orderElem + 1, 0.0);

        /*1. Calculate volume integral term
         * Chi excute khi order > 0
        */
        if (mathVar::orderElem>0)
        {
            process::auxEq::BR1::calcVolumeIntegralTerms(element, rhoVolIntOx, rhouVolIntOx, rhovVolIntOx, rhoEVolIntOx, rhoVolIntOy, rhouVolIntOy, rhovVolIntOy, rhoEVolIntOy);
        }

        /*2. Calculate surface integral term*/
        process::auxEq::BR1::calcSurfaceIntegralTerms(element, rhoSurfIntOx, rhouSurfIntOx, rhovSurfIntOx, rhoESurfIntOx, rhoSurfIntOy, rhouSurfIntOy, rhovSurfIntOy, rhoESurfIntOy);

        for (int order = 0; order <= mathVar::orderElem; order++)
        {
            rhoRHSOx[order] = -rhoVolIntOx[order] + rhoSurfIntOx[order];
            rhoRHSOy[order] = -rhoVolIntOy[order] + rhoSurfIntOy[order];

            rhouRHSOx[order] = -rhouVolIntOx[order] + rhouSurfIntOx[order];
            rhouRHSOy[order] = -rhouVolIntOy[order] + rhouSurfIntOy[order];

            rhovRHSOx[order] = -rhovVolIntOx[order] + rhovSurfIntOx[order];
            rhovRHSOy[order] = -rhovVolIntOy[order] + rhovSurfIntOy[order];

            rhoERHSOx[order] = -rhoEVolIntOx[order] + rhoESurfIntOx[order];
            rhoERHSOy[order] = -rhoEVolIntOy[order] + rhoESurfIntOy[order];
        }
    }

    void calcVolumeIntegralTerms(int element, std::vector<double> &rhoVolIntX, std::vector<double> &rhouVolIntX, std::vector<double> &rhovVolIntX, std::vector<double> &rhoEVolIntX, std::vector<double> &rhoVolIntY, std::vector<double> &rhouVolIntY, std::vector<double> &rhovVolIntY, std::vector<double> &rhoEVolIntY)
    {
        std::vector<std::vector<double>> rhoGsVol(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
            rhouGsVol(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
            rhovGsVol(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
            rhoEGsVol(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));

        //rho -------------------------------------------------------------------------------------------
        process::auxEq::BR1::getGaussMatrixOfConserVar(rhoGsVol, element, 1);
        //rhou ------------------------------------------------------------------------------------------
        process::auxEq::BR1::getGaussMatrixOfConserVar(rhouGsVol, element, 2);
        //rhov ------------------------------------------------------------------------------------------
        process::auxEq::BR1::getGaussMatrixOfConserVar(rhovGsVol, element, 3);
        //rhoE ------------------------------------------------------------------------------------------
        process::auxEq::BR1::getGaussMatrixOfConserVar(rhoEGsVol, element, 4);

        //Bien phu la mu*divU nen can nhan mu truoc khi tinh tich phan
        double mu(0.0);
        for (int na = 0; na <= mathVar::nGauss2D; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
            {
                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                mu=math::CalcVisCoef(volumeFields::T[element][nanb]);
                rhoGsVol[na][nb]=rhoGsVol[na][nb]*mu;
                rhouGsVol[na][nb]=rhouGsVol[na][nb]*mu;
                rhovGsVol[na][nb]=rhovGsVol[na][nb]*mu;
                rhoEGsVol[na][nb]=rhoEGsVol[na][nb]*mu;
            }
        }

        //Tinh tich phan
        for (int order = 1; order <= mathVar::orderElem; order++)
        {
            rhoVolIntX[order] = process::volumeInte(element, rhoGsVol, order, 1);
            rhoVolIntY[order] = process::volumeInte(element, rhoGsVol, order, 2);

            rhouVolIntX[order] = process::volumeInte(element, rhouGsVol, order, 1);
            rhovVolIntX[order] = process::volumeInte(element, rhovGsVol, order, 1);
            rhoEVolIntX[order] = process::volumeInte(element, rhoEGsVol, order, 1);

            rhouVolIntY[order] = process::volumeInte(element, rhouGsVol, order, 2);
            rhovVolIntY[order] = process::volumeInte(element, rhovGsVol, order, 2);
            rhoEVolIntY[order] = process::volumeInte(element, rhoEGsVol, order, 2);
        }
        rhoVolIntX[0] = 0;
        rhoVolIntY[0] = 0;

        rhouVolIntX[0] = 0;
        rhovVolIntX[0] = 0;
        rhoEVolIntX[0] = 0;

        rhouVolIntY[0] = 0;
        rhovVolIntY[0] = 0;
        rhoEVolIntY[0] = 0;
    }

    void calcSurfaceIntegralTerms(int element, std::vector<double> &rhoSurfIntX, std::vector<double> &rhouSurfIntX, std::vector<double> &rhovSurfIntX, std::vector<double> &rhoESurfIntX, std::vector<double> &rhoSurfIntY, std::vector<double> &rhouSurfIntY, std::vector<double> &rhovSurfIntY, std::vector<double> &rhoESurfIntY)
    {
        /*User's guide:
        Input array rhoRHSTerm, rhouRHSTerm, rhovRHSTerm, rhoERHSTerm have following form:
        - number of row: orderElem + 1*/
        int elemType(auxUlti::checkType(element)), edgeName(0);
        std::vector<std::vector<double>> rhoFluxX(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
            rhouFluxX(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
            rhovFluxX(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
            rhoEFluxX(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
            rhoFluxY(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
            rhouFluxY(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
            rhovFluxY(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
            rhoEFluxY(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0));

        std::vector<double> rhoFluxXTemp(mathVar::nGauss1D + 1, 0.0),
            rhouFluxXTemp(mathVar::nGauss1D + 1, 0.0),
            rhovFluxXTemp(mathVar::nGauss1D + 1, 0.0),
            rhoEFluxXTemp(mathVar::nGauss1D + 1, 0.0),
            rhoFluxYTemp(mathVar::nGauss1D + 1, 0.0),
            rhouFluxYTemp(mathVar::nGauss1D + 1, 0.0),
            rhovFluxYTemp(mathVar::nGauss1D + 1, 0.0),
            rhoEFluxYTemp(mathVar::nGauss1D + 1, 0.0);

        /*1. Calculate flux of conservative variables at all Gauss points on all faces of element*/
        process::auxEq::BR1::getGaussVectorOfConserVar(element, rhoFluxX, rhouFluxX, rhovFluxX, rhoEFluxX, rhoFluxY, rhouFluxY, rhovFluxY, rhoEFluxY);

        /*2. Calculates surface integrals of all conservative variables at all order*/
        for (int nface = 0; nface < elemType; nface++)
        {
            edgeName = meshVar::inedel[element][nface];
            for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
            {
                rhoFluxXTemp[nG] = rhoFluxX[nG][nface];
                rhoFluxYTemp[nG] = rhoFluxY[nG][nface];

                rhouFluxXTemp[nG] = rhouFluxX[nG][nface];
                rhovFluxXTemp[nG] = rhovFluxX[nG][nface];
                rhoEFluxXTemp[nG] = rhoEFluxX[nG][nface];

                rhouFluxYTemp[nG] = rhouFluxY[nG][nface];
                rhovFluxYTemp[nG] = rhovFluxY[nG][nface];
                rhoEFluxYTemp[nG] = rhoEFluxY[nG][nface];
            }

            for (int order = 0; order <= mathVar::orderElem; order++)
            {
                rhoSurfIntX[order] += process::surfaceInte(element, edgeName, rhoFluxXTemp, order);
                rhoSurfIntY[order] += process::surfaceInte(element, edgeName, rhoFluxYTemp, order);

                rhouSurfIntX[order] += process::surfaceInte(element, edgeName, rhouFluxXTemp, order);
                rhovSurfIntX[order] += process::surfaceInte(element, edgeName, rhovFluxXTemp, order);
                rhoESurfIntX[order] += process::surfaceInte(element, edgeName, rhoEFluxXTemp, order);

                rhouSurfIntY[order] += process::surfaceInte(element, edgeName, rhouFluxYTemp, order);
                rhovSurfIntY[order] += process::surfaceInte(element, edgeName, rhovFluxYTemp, order);
                rhoESurfIntY[order] += process::surfaceInte(element, edgeName, rhoEFluxYTemp, order);
            }
        }
    }

    void getGaussMatrixOfConserVar(std::vector<std::vector<double>> &UGauss, int element, int valType)
    {
        //rhoE used for solving auxilary equation is advective rhoE!!
        //std::vector<std::vector<double>> GaussMatrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
        for (int na = 0; na <= mathVar::nGauss2D; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
            {
                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                switch (valType)
                {
                case 1:
                {
                    UGauss[na][nb] = volumeFields::rhoVolGauss[element][nanb];
                }
                break;
                case 2:
                {
                    UGauss[na][nb] = volumeFields::rhouVolGauss[element][nanb];
                }
                break;
                case 3:
                {
                    UGauss[na][nb] = volumeFields::rhovVolGauss[element][nanb];
                }
                break;
                case 4:
                {
                    UGauss[na][nb] = volumeFields::rhoEVolGauss[element][nanb];
                }
                break;
                default:
                    break;
                }
            }
        }
    }

    std::vector<std::vector<double>> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nx, double ny)
    {
        //NOTE: Ham nay dung khi dai phuong trinh phu, vi ham nay tra ve gia tri +/- cua mu*U
        std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0)); //columns 0, 1 are Ox, Oy values
        double rhoM(0.0), rhoP(0.0), rhouM(0.0), rhouP(0.0), rhovM(0.0), rhovP(0.0), TP(0.0), TM(0.0), rhoEM(0.0), rhoEP(0.0), muM, muP;
        std::tie(rhoP, rhoM) = auxUlti::getUAtInterfaces(edge, element, nG, 1);
        std::tie(rhouP, rhouM) = auxUlti::getUAtInterfaces(edge, element, nG, 2);
        std::tie(rhovP, rhovM) = auxUlti::getUAtInterfaces(edge, element, nG, 3);

        std::tie(TP, TM)=auxUlti::getTAtInterfaces(edge,element,nG);
        muP=math::CalcVisCoef(TP);
        muM=math::CalcVisCoef(TM);

        std::tie(rhoEP, rhoEM) = auxUlti::getUAtInterfaces(edge, element, nG, 4);

        //Tinh Flux
        gaussVector[0][0] = math::numericalFluxes::auxFlux(rhoM*muM, rhoP*muP, nx);
        gaussVector[0][1] = math::numericalFluxes::auxFlux(rhoM*muM, rhoP*muP, ny);
        gaussVector[1][0] = math::numericalFluxes::auxFlux(rhouM*muM, rhouP*muP, nx);
        gaussVector[1][1] = math::numericalFluxes::auxFlux(rhouM*muM, rhouP*muP, ny);
        gaussVector[2][0] = math::numericalFluxes::auxFlux(rhovM*muM, rhovP*muP, nx);
        gaussVector[2][1] = math::numericalFluxes::auxFlux(rhovM*muM, rhovP*muP, ny);
        gaussVector[3][0] = math::numericalFluxes::auxFlux(rhoEM*muM, rhoEP*muP, nx);
        gaussVector[3][1] = math::numericalFluxes::auxFlux(rhoEM*muM, rhoEP*muP, ny);

        return gaussVector;
    }

    void getGaussVectorOfConserVar(int element, std::vector<std::vector<double>> &rhoFluxX, std::vector<std::vector<double>> &rhouFluxX, std::vector<std::vector<double>> &rhovFluxX, std::vector<std::vector<double>> &rhoEFluxX, std::vector<std::vector<double>> &rhoFluxY, std::vector<std::vector<double>> &rhouFluxY, std::vector<std::vector<double>> &rhovFluxY, std::vector<std::vector<double>> &rhoEFluxY)
    {
        /*User's guide:
        Input array rhoFlux, rhouFlux, rhovFlux, rhoEFlux have following form:
        - number of row: nGauss + 1*/

        int elemType(auxUlti::checkType(element)), edgeName(0);
        int faceBcType(0);
        std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0));

        for (int nface = 0; nface < elemType; nface++)
        {
            edgeName = meshVar::inedel[element][nface];
            faceBcType = auxUlti::getBCType(edgeName);
            double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2));
            if (faceBcType == 0)  //internal edge
            {
                for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
                {
                    gaussVector = process::auxEq::BR1::getVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss, nx, ny);
                    rhoFluxX[nGauss][nface] = gaussVector[0][0];
                    rhouFluxX[nGauss][nface] = gaussVector[1][0];
                    rhovFluxX[nGauss][nface] = gaussVector[2][0];
                    rhoEFluxX[nGauss][nface] = gaussVector[3][0];

                    rhoFluxY[nGauss][nface] = gaussVector[0][1];
                    rhouFluxY[nGauss][nface] = gaussVector[1][1];
                    rhovFluxY[nGauss][nface] = gaussVector[2][1];
                    rhoEFluxY[nGauss][nface] = gaussVector[3][1];
                }
            }
            else  //boundary edge
            {
                for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
                {
                    gaussVector = auxEqBCsImplement(element, edgeName, nGauss);
                    rhoFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[0][0],gaussVector[0][1],nx);
                    rhouFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[1][0],gaussVector[1][1],nx);
                    rhovFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[2][0],gaussVector[2][1],nx);
                    rhoEFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[3][0],gaussVector[3][1],nx);

                    rhoFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[0][0],gaussVector[0][1],ny);
                    rhouFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[1][0],gaussVector[1][1],ny);
                    rhovFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[2][0],gaussVector[2][1],ny);
                    rhoEFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[3][0],gaussVector[3][1],ny);
                }
            }
        }
    }
    }

	}//end namespace auxEq

	namespace NSFEq
	{
        void calcLocalInviscidFlux(std::vector<double> &UG, std::vector<double> &FLocalX, std::vector<double> &FLocalY, std::vector<double> &UL, std::vector<double> &n, double T)
        {
            /* Ham tinh inviscous term tu U,
             * - Xoay U ve he toa do local truoc khi tinh
             * - Gia tri xuat ra la gia tri F_invs da xoay ve he local
             * - Chu y: neu flux type la Lax Friedrichs thi k can xoay
             *
             * Cac bien tra ve (thong qua pass by reference):
             * - FlocalX: inviscid flux theo phuong Ox local, hay la phuong normal
             * - FlocalY: inviscid flux theo phuong Oy local, hay la phuong tangential
             * - UL: conservative variable da duoc xoay ve he local
            */

            double nx(n[0]), ny(n[1]);
            double rho(UG[0]), rhouG(UG[1]), rhovG(UG[2]), rhoE(UG[3]), rhouG_M, rhovG_M; //G = global, M = mass diffusion
            double rhouL, //L = local, rhouL tuong duong normal direction
                    rhovL, //rhovL tuong duong tangential direction
                    rhouL_M,
                    rhovL_M;

            rhouG_M=rhouG;
            rhovG_M=rhovG;

            //Xoay ve he local
            if (DGSchemes::fluxControl::LxF)
            {
                rhouL=rhouG;
                rhovL=rhovG;
                rhouL_M=rhouG_M;
                rhovL_M=rhovG_M;
            }
            else
            {
                std::tie(rhouL,rhovL)=math::rotateToLocalCoordinateSystem(rhouG,rhovG,nx,ny);
                std::tie(rhouL_M,rhovL_M)=math::rotateToLocalCoordinateSystem(rhouG_M,rhovG_M,nx,ny);
            }

            //Cap nhat vector UL
            /* Neu mass diffusion ON, UL[1] va UL[2] la total momentum.
             * Neu mass diffusion OFF, UL[1] va UL[2] la convective momentum
            */
            UL[0]=rho;
            UL[1]=rhouL;
            UL[2]=rhovL;
            UL[3]=rhoE;

            //Tinh inviscid flux
            double
            uL = (rhouL / rho),
            vL = (rhovL / rho),
            uL_M = (rhouL_M / rho),
            vL_M = (rhovL_M / rho),
            totalE = (rhoE / rho),

            /*calculate P*/
            p = math::CalcP(T, rho);

            /*calculate inviscid terms*/
            std::tie(FLocalX[0], FLocalX[1], FLocalX[2], FLocalX[3]) = math::inviscidTerms::calcInvisTermsFromPriVars(rho, uL, uL_M, vL, vL_M, totalE, p, 1);
            std::tie(FLocalY[0], FLocalY[1], FLocalY[2], FLocalY[3]) = math::inviscidTerms::calcInvisTermsFromPriVars(rho, uL, uL_M, vL, vL_M, totalE, p, 2);
        }

        void calcLocalViscousFlux(std::vector<double> &UG, std::vector<double> &dUX, std::vector<double> &dUY, std::vector<double> &FLocalX, std::vector<double> &FLocalY, std::vector<double> &n)
        {
            /* Tai version hien tai, diffusive flux va auxilarity flux deu su dung central flux, vi vay khong can xoay ve he local,
             * tuy nhien trong code nay van khai bao cac bien local de co the de dang cap nhat cac flux type khac cho diffusive flux ve sau.
            */

            std::vector<double>UL(4,0.0), //L = local
                    dULocalX(4,0.0), //norm tuong duong Ox
                    dULocalY(4,0.0); //tang tuong duong Oy
            //double nx(n[0]), ny(n[1]);

            //Xoay ve he local
            //UL[0]=UG[0];
            //std::tie(UL[1],UL[2])=math::rotateToLocalCoordinateSystem(UG[1],UG[2],nx,ny);
            //UL[3]=UG[3];
            /*
            for (int i=0; i<4; i++)
            {
                std::tie(dULocalX[i],dULocalY[i])=math::rotateToLocalCoordinateSystem(dUX[i],dUY[i],nx,ny);
            }
            */

            UL=UG;
            dULocalX=dUX;
            dULocalY=dUY;

            /*StressHeat matrix has form:
            [tauXx		tauXy		Qx]
            [tauYx		tauYy		Qy]
            */
            std::vector<std::vector<double>> StressHeat(2, std::vector<double>(3, 0.0));

            StressHeat = math::viscousTerms::calcStressTensorAndHeatFlux(UL, dULocalX, dULocalY);
            std::tie(FLocalX[0], FLocalX[1], FLocalX[2], FLocalX[3]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeat, UL[1]/UL[0], UL[2]/UL[0], 1);
            std::tie(FLocalY[0], FLocalY[1], FLocalY[2], FLocalY[3]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeat, UL[1]/UL[0], UL[2]/UL[0], 2);

            //DURST MODEL ------------------------------------------------------------------------------
            //Correct Diffusive Term of NSF Eqns
            //Theo mo hinh Durst
            if (extNSF_Durst::enable)
            {
                std::vector<std::vector<double>> diffTerms(2, std::vector<double>(4, 0.0));

                //Chi khi nao mode diffusionAtWall la No (khong cho phep diffusion at wall)
                //Va flag needToRemoveDiffTerm dang True thi moi remove Diff Terms
                if (!extNSF_Durst::diffusionAtWall && extNSF_Durst::needToRemoveDiffTerm)
                {
                    //Do nothing
                }
                else if (
                         //Chi khi nao mode diffusionAtWall la Yes (cho phep diffusion at wall)
                         //Va flag dropNormSelfDiffTerm dang True thi moi remove normal Diff Terms (bang cach dung ham bcForExtNSF_Durst::dropNormSelfDiffTerm)
                         (extNSF_Durst::diffusionAtWall && extNSF_Durst::dropNormSelfDiffTerm)

                         //Hoac BC dang tinh flux la symmetry
                         || (extNSF_Durst::isSymmetry)
                         )
                {
                    bcForExtNSF_Durst::correctViscousTerms(diffTerms,UL,dULocalX,dULocalY,n);
                }
                //Cac truong hop khac: volume term va internal edge
                else
                {
                    extNSF_Durst::correctViscousTerms(diffTerms,UL,dULocalX,dULocalY);
                }

                //Add to classical diffusive terms
                for (int i=0; i<4; i++)
                {
                    FLocalX[i] = FLocalX[i] + diffTerms[0][i];
                    FLocalY[i] = FLocalY[i] + diffTerms[1][i];
                }
            }
            //DURST MODEL ------------------------------------------------------------------------------
        }

        void calcSurfaceFlux()
        {
            //Tinh flux qua cac internal edge.
            //mu included
            int masterCell(-1), slaveCell(-1), bcGrp(0);
            double TMaster(0.0), TSlave(0.0);

            //Global coordinates
            std::vector<double> UMaster(4, 0.0), dUXMaster(4, 0.0), dUYMaster(4, 0.0),
                USlave(4, 0.0), dUXSlave(4, 0.0), dUYSlave(4, 0.0), vectorn(2, 0.0);

            //Local coordinates
            std::vector<double> ULMaster(4, 0.0), ULSlave(4, 0.0);

            std::vector<double>
                    //Local fluxes
                    invisFluxLocalXMaster(4,0.0),
                    invisFluxLocalYMaster(4,0.0),
                    visFluxLocalXMaster(4,0.0),
                    visFluxLocalYMaster(4,0.0),

                    invisFluxLocalXSlave(4,0.0),
                    invisFluxLocalYSlave(4,0.0),
                    visFluxLocalXSlave(4,0.0),
                    visFluxLocalYSlave(4,0.0),

                    //vector lay gia tri flux o he toa do global
                    invisFluxGlobal(4,0.0),
                    visFluxGlobal(4,0.0);

            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
            {
                bcGrp = auxUlti::getBCType(iedge);
                if (bcGrp == 0)
                {
                    //DURST MODEL ------------------------------------------------------------------------------
                    extNSF_Durst::needToRemoveDiffTerm=false;
                    //DURST MODEL ------------------------------------------------------------------------------

                    std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    vectorn[0] = auxUlti::getNormVectorComp(masterCell, iedge, 1);
                    vectorn[1] = auxUlti::getNormVectorComp(masterCell, iedge, 2);

                    if (DGSchemes::fluxControl::LxF)
                    {
                        //Tinh he so C cua LxF flux
                        math::numericalFluxes::findMaxLxFConstantOnEdge(iedge,masterCell);
                    }

                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        std::tie(TMaster,TSlave)=auxUlti::getTAtInterfaces(iedge,masterCell,nG);
                        for (int i = 0; i < 4; i++)
                        {
                            //std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
                            std::tie(UMaster[i], USlave[i]) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, i + 1);
                            //d(U)/dx, d(U)/dy
                            if (flowProperties::viscous)
                            {
                                std::tie(dUXMaster[i], dUXSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i+1, 1);
                                std::tie(dUYMaster[i], dUYSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i+1, 2);
                            }
                        }

                        if (flowProperties::viscous)
                        {
                            //Save derivaties to surfaceFields
                            surfaceFields::dRhoX[iedge][nG]=dUXMaster[0];
                            surfaceFields::dRhouX[iedge][nG]=dUXMaster[1];
                            surfaceFields::dRhovX[iedge][nG]=dUXMaster[2];
                            surfaceFields::dRhoEX[iedge][nG]=dUXMaster[3];

                            surfaceFields::dRhoX[iedge][mathVar::nGauss1D+nG+1]=dUXSlave[0];
                            surfaceFields::dRhouX[iedge][mathVar::nGauss1D+nG+1]=dUXSlave[1];
                            surfaceFields::dRhovX[iedge][mathVar::nGauss1D+nG+1]=dUXSlave[2];
                            surfaceFields::dRhoEX[iedge][mathVar::nGauss1D+nG+1]=dUXSlave[3];

                            surfaceFields::dRhoY[iedge][nG]=dUYMaster[0];
                            surfaceFields::dRhouY[iedge][nG]=dUYMaster[1];
                            surfaceFields::dRhovY[iedge][nG]=dUYMaster[2];
                            surfaceFields::dRhoEY[iedge][nG]=dUYMaster[3];

                            surfaceFields::dRhoY[iedge][mathVar::nGauss1D+nG+1]=dUYSlave[0];
                            surfaceFields::dRhouY[iedge][mathVar::nGauss1D+nG+1]=dUYSlave[1];
                            surfaceFields::dRhovY[iedge][mathVar::nGauss1D+nG+1]=dUYSlave[2];
                            surfaceFields::dRhoEY[iedge][mathVar::nGauss1D+nG+1]=dUYSlave[3];
                        }

                        /* Voi ham calcLocalInviscidFlux, neu flux type la LxF hoac central, cac bien
                         * invisFluxLocalXMaster,.., ULMaster,.. deu luu gia tri tai he toa do global.
                         * Neu flux type la Roe hoac HLL, cac bien nay luu gia tri tai he local.
                        */
                        process::NSFEq::calcLocalInviscidFlux(UMaster,invisFluxLocalXMaster,invisFluxLocalYMaster,ULMaster,vectorn,TMaster);
                        process::NSFEq::calcLocalInviscidFlux(USlave,invisFluxLocalXSlave,invisFluxLocalYSlave,ULSlave,vectorn,TSlave);

                        /* Voi ham calcLocalViscousFlux, vi flux type luon la central, cac bien
                         * invisFluxLocalXMaster,.., ULMaster,.. deu luu gia tri tai he toa do global.
                        */
                        if (flowProperties::viscous)
                        {
                            process::NSFEq::calcLocalViscousFlux(UMaster,dUXMaster,dUYMaster,visFluxLocalXMaster,visFluxLocalYMaster,vectorn);
                            process::NSFEq::calcLocalViscousFlux(USlave,dUXSlave,dUYSlave,visFluxLocalXSlave,visFluxLocalYSlave,vectorn);
                        }

                        //Tinh flux
                        //Vi dang tinh cho cell master cua edge nen phia trong (phia +) la cua cell master
                        math::numericalFluxes::calConvectiveFluxes(iedge, invisFluxGlobal,
                                                               invisFluxLocalXMaster, invisFluxLocalXSlave,
                                                               invisFluxLocalYMaster, invisFluxLocalYSlave,
                                                               ULMaster, ULSlave,
                                                               TMaster, TSlave,
                                                               vectorn);

                        //Viscous flux tinh bang central flux
                        if (flowProperties::viscous)
                        {
                            math::numericalFluxes::centralFlux(visFluxGlobal,
                                                               visFluxLocalXMaster, visFluxLocalXSlave,
                                                               visFluxLocalYMaster, visFluxLocalYSlave,
                                                               vectorn);
                        }

                        //Luu flux vao vector
                        /* Flux tren 1 mat co tinh bao toan (conservative) nen fluxMaster = -fluxSlave
                        */
                        surfaceFields::invis_rho[iedge][nG]=invisFluxGlobal[0];
                        surfaceFields::invis_rhou[iedge][nG]=invisFluxGlobal[1];
                        surfaceFields::invis_rhov[iedge][nG]=invisFluxGlobal[2];
                        surfaceFields::invis_rhoE[iedge][nG]=invisFluxGlobal[3];

                        surfaceFields::invis_rho[iedge][nG + mathVar::nGauss1D + 1]=-invisFluxGlobal[0];
                        surfaceFields::invis_rhou[iedge][nG + mathVar::nGauss1D + 1]=-invisFluxGlobal[1];
                        surfaceFields::invis_rhov[iedge][nG + mathVar::nGauss1D + 1]=-invisFluxGlobal[2];
                        surfaceFields::invis_rhoE[iedge][nG + mathVar::nGauss1D + 1]=-invisFluxGlobal[3];

                        if (flowProperties::viscous)
                        {
                            surfaceFields::Vis_rho[iedge][nG]=visFluxGlobal[0];
                            surfaceFields::Vis_rhou[iedge][nG]=visFluxGlobal[1];
                            surfaceFields::Vis_rhov[iedge][nG]=visFluxGlobal[2];
                            surfaceFields::Vis_rhoE[iedge][nG]=visFluxGlobal[3];

                            surfaceFields::Vis_rho[iedge][nG + mathVar::nGauss1D + 1]=-visFluxGlobal[0];
                            surfaceFields::Vis_rhou[iedge][nG + mathVar::nGauss1D + 1]=-visFluxGlobal[1];
                            surfaceFields::Vis_rhov[iedge][nG + mathVar::nGauss1D + 1]=-visFluxGlobal[2];
                            surfaceFields::Vis_rhoE[iedge][nG + mathVar::nGauss1D + 1]=-visFluxGlobal[3];
                        }
                    }
                }
                /*
                else if (bcGrp == 4)
                {
                    if (flowProperties::viscous)
                    {
                        int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(iedge));
                        std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);

                        for (int nG = 0; nG <= mathVar::nGauss; nG++)
                        {
                            SurfaceBCFields::GaussDRhoX[loc][nG] = BR1Vars::rhoX[masterCell][0];
                            SurfaceBCFields::GaussDRhouX[loc][nG] = BR1Vars::rhouX[masterCell][0];
                            SurfaceBCFields::GaussDRhovX[loc][nG] = BR1Vars::rhovX[masterCell][0];
                            SurfaceBCFields::GaussDRhoEX[loc][nG] = BR1Vars::rhoEX[masterCell][0];

                            SurfaceBCFields::GaussDRhoY[loc][nG] = BR1Vars::rhoY[masterCell][0];
                            SurfaceBCFields::GaussDRhouY[loc][nG] = BR1Vars::rhouY[masterCell][0];
                            SurfaceBCFields::GaussDRhovY[loc][nG] = BR1Vars::rhovY[masterCell][0];
                            SurfaceBCFields::GaussDRhoEY[loc][nG] = BR1Vars::rhoEY[masterCell][0];
                        }
                    }
                }*/
                else
                {
                    if (flowProperties::viscous)
                    {
                        //int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(iedge));
                        std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);

                        for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                        {
                            surfaceFields::dRhoX[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,1,1);
                            surfaceFields::dRhouX[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,2,1);
                            surfaceFields::dRhovX[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,3,1);
                            surfaceFields::dRhoEX[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,4,1);

                            surfaceFields::dRhoY[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,1,2);
                            surfaceFields::dRhouY[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,2,2);
                            surfaceFields::dRhovY[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,3,2);
                            surfaceFields::dRhoEY[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,4,2);
                        }
                    }
                }
            }

            if (systemVar::parallelMode && flowProperties::viscous)
            {
                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhoX, surfaceFields::dRhoX);
                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhouX, surfaceFields::dRhouX);
                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhovX, surfaceFields::dRhovX);
                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhoEX, surfaceFields::dRhoEX);

                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhoY, surfaceFields::dRhoY);
                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhouY, surfaceFields::dRhouY);
                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhovY, surfaceFields::dRhovY);
                parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhoEY, surfaceFields::dRhoEY);
            }
        }

        void solveNSFEquation_TVDRK(int RKOrder)
		{
			std::vector<double> rhoError(meshVar::nelem2D, 1.0),
				rhouError(meshVar::nelem2D, 1.0),
				rhovError(meshVar::nelem2D, 1.0),
				rhoEError(meshVar::nelem2D, 1.0);

			//std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
			std::vector<double>
				RHSTerm1(mathVar::orderElem + 1, 0.0),
				RHSTerm2(mathVar::orderElem + 1, 0.0),
				RHSTerm3(mathVar::orderElem + 1, 0.0),
				RHSTerm4(mathVar::orderElem + 1, 0.0),

				ddtRhoVector(mathVar::orderElem + 1, 0.0),
				ddtRhouVector(mathVar::orderElem + 1, 0.0),
				ddtRhovVector(mathVar::orderElem + 1, 0.0),
				ddtRhoEVector(mathVar::orderElem + 1, 0.0),

				rhoVectorN(mathVar::orderElem + 1, 0.0),
				rhouVectorN(mathVar::orderElem + 1, 0.0),
				rhovVectorN(mathVar::orderElem + 1, 0.0),
				rhoEVectorN(mathVar::orderElem + 1, 0.0),

				UnVector(mathVar::orderElem + 1, 0.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				//2) Calculate Right hand side terms
				process::NSFEq::CalcRHSTerm(nelement, RHSTerm1, RHSTerm2, RHSTerm3, RHSTerm4);

                //3) Solve for time derivartives of conservative variables
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					ddtRhoVector[iorder] = RHSTerm1[iorder] / stiffMatrixCoeffs[nelement][iorder];
					ddtRhouVector[iorder] = RHSTerm2[iorder] / stiffMatrixCoeffs[nelement][iorder];
					ddtRhovVector[iorder] = RHSTerm3[iorder] / stiffMatrixCoeffs[nelement][iorder];
					ddtRhoEVector[iorder] = RHSTerm4[iorder] / stiffMatrixCoeffs[nelement][iorder];

					switch (RKOrder)
					{
					case 1:
					{
						rhoResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhoVector[iorder];
						rhouResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhouVector[iorder];
						rhovResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhovVector[iorder];
						rhoEResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhoEVector[iorder];
					}
					break;
					case 2:
					{
						rhoResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhoVector[iorder];
						rhouResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhouVector[iorder];
						rhovResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhovVector[iorder];
						rhoEResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhoEVector[iorder];
					}
					break;
					case 3:
					{
						rhoResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhoVector[iorder];
						rhouResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhouVector[iorder];
						rhovResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhovVector[iorder];
						rhoEResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhoEVector[iorder];
					}
					break;
					default:
						break;
					}
				}


				//4) Solve time marching
				//rho
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rho[nelement][order];
				}
				rhoVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhoVector, UnVector, RKOrder, 1);
				//rhou
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhou[nelement][order];
				}
				rhouVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhouVector, UnVector, RKOrder, 2);
				//rhov
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhov[nelement][order];
				}
				rhovVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhovVector, UnVector, RKOrder, 3);
				//rhoE
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhoE[nelement][order];
				}
				rhoEVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhoEVector, UnVector, RKOrder, 4);

                //5) Save results to conservative variables arrays
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoN[nelement][order] = rhoVectorN[order];
					rhouN[nelement][order] = rhouVectorN[order];
					rhovN[nelement][order] = rhovVectorN[order];
					rhoEN[nelement][order] = rhoEVectorN[order];
				}
			}
		}

        void solveNSFEquation_Euler()
        {
            std::vector<double> rhoError(meshVar::nelem2D, 1.0),
                rhouError(meshVar::nelem2D, 1.0),
                rhovError(meshVar::nelem2D, 1.0),
                rhoEError(meshVar::nelem2D, 1.0);

            //std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
            std::vector<double>
                RHSTerm1(mathVar::orderElem + 1, 0.0),
                RHSTerm2(mathVar::orderElem + 1, 0.0),
                RHSTerm3(mathVar::orderElem + 1, 0.0),
                RHSTerm4(mathVar::orderElem + 1, 0.0),

                ddtRhoVector(mathVar::orderElem + 1, 0.0),
                ddtRhouVector(mathVar::orderElem + 1, 0.0),
                ddtRhovVector(mathVar::orderElem + 1, 0.0),
                ddtRhoEVector(mathVar::orderElem + 1, 0.0),

                rhoVectorN(mathVar::orderElem + 1, 0.0),
                rhouVectorN(mathVar::orderElem + 1, 0.0),
                rhovVectorN(mathVar::orderElem + 1, 0.0),
                rhoEVectorN(mathVar::orderElem + 1, 0.0),

                UnVector(mathVar::orderElem + 1, 0.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
            {
                //2) Calculate Right hand side terms
                process::NSFEq::CalcRHSTerm(nelement, RHSTerm1, RHSTerm2, RHSTerm3, RHSTerm4);

                //3) Solve for time derivartives of conservative variables
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    ddtRhoVector[iorder] = RHSTerm1[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    ddtRhouVector[iorder] = RHSTerm2[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    ddtRhovVector[iorder] = RHSTerm3[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    ddtRhoEVector[iorder] = RHSTerm4[iorder] / stiffMatrixCoeffs[nelement][iorder];

                    rhoResArr[nelement][iorder] = ddtRhoVector[iorder];
                    rhouResArr[nelement][iorder] = ddtRhouVector[iorder];
                    rhovResArr[nelement][iorder] = ddtRhovVector[iorder];
                    rhoEResArr[nelement][iorder] = ddtRhoEVector[iorder];
                }

                //4) Solve time marching
                //rho
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    UnVector[order] = rho[nelement][order];
                }
                rhoVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhoVector, UnVector, 1, 1);
                //rhou
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    UnVector[order] = rhou[nelement][order];
                }
                rhouVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhouVector, UnVector, 1, 2);
                //rhov
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    UnVector[order] = rhov[nelement][order];
                }
                rhovVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhovVector, UnVector, 1, 3);
                //rhoE
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    UnVector[order] = rhoE[nelement][order];
                }
                rhoEVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhoEVector, UnVector, 1, 4);

                //5) Save results to conservative variables arrays
                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    rhoN[nelement][order] = rhoVectorN[order];
                    rhouN[nelement][order] = rhouVectorN[order];
                    rhovN[nelement][order] = rhovVectorN[order];
                    rhoEN[nelement][order] = rhoEVectorN[order];
                }
            }
        }

		void updateVariables()
		{
            rho = rhoN;
            rhou = rhouN;
            rhov = rhovN;
            rhoE = rhoEN;
		}

		/*Function calculates right hand side terms of all conservative variables at ONLY one order*/
		void CalcRHSTerm(int element, std::vector<double> &term1RHS, std::vector<double> &term2RHS, std::vector<double> &term3RHS, std::vector<double> &term4RHS)
		{
			std::vector<double>
				VolIntTerm1(mathVar::orderElem + 1, 0.0),
				VolIntTerm2(mathVar::orderElem + 1, 0.0),
				VolIntTerm3(mathVar::orderElem + 1, 0.0),
				VolIntTerm4(mathVar::orderElem + 1, 0.0),

				SurfIntTerm1(mathVar::orderElem + 1, 0.0),
				SurfIntTerm2(mathVar::orderElem + 1, 0.0),
				SurfIntTerm3(mathVar::orderElem + 1, 0.0),
				SurfIntTerm4(mathVar::orderElem + 1, 0.0);

			/*Volume integral term===========================================================================*/
            if (mathVar::orderElem>0)
                process::NSFEq::calcVolumeIntegralTerms(element, VolIntTerm1, VolIntTerm2, VolIntTerm3, VolIntTerm4);
			
			/*Surface integral term===========================================================================*/
			process::NSFEq::calcSurfaceIntegralTerms(element, SurfIntTerm1, SurfIntTerm2, SurfIntTerm3, SurfIntTerm4);
			
            for (int order = 0; order <= mathVar::orderElem; order++)
			{
				term1RHS[order] = VolIntTerm1[order] - SurfIntTerm1[order];
				term2RHS[order] = VolIntTerm2[order] - SurfIntTerm2[order];
				term3RHS[order] = VolIntTerm3[order] - SurfIntTerm3[order];
				term4RHS[order] = VolIntTerm4[order] - SurfIntTerm4[order];
			}
		}

		/*Function calculates Inviscid terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussInviscidTerm(int element, int na, int nb)
		{
			/*InviscidTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
            int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
			std::vector<std::vector<double>> InviscidTerm(4, std::vector<double>(2, 0.0));
			double
                rhoVal(volumeFields::rhoVolGauss[element][nanb]),
                rhouVal(volumeFields::rhouVolGauss[element][nanb]),
                rhovVal(volumeFields::rhovVolGauss[element][nanb]),
                rhoEVal(volumeFields::rhoEVolGauss[element][nanb]),
				uVal(0.0),
				vVal(0.0),
				pVal(0.0),
                umVal(0.0),
                vmVal(0.0),
                totalE(0.0),
                a(0.0),b(0.0);//muVal(0.0);

            std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
			uVal = rhouVal / rhoVal;
			vVal = rhovVal / rhoVal;
            totalE = rhoEVal / rhoVal;
            pVal = math::CalcP(volumeFields::T[element][nanb], rhoVal);
            //muVal=math::CalcVisCoef(volumeFields::T[element][nanb]);

            umVal=uVal;
            vmVal=vVal;

            /*1. Ox direction*/
            std::tie(InviscidTerm[0][0], InviscidTerm[1][0], InviscidTerm[2][0], InviscidTerm[3][0]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, umVal, vVal, vmVal, totalE, pVal, 1);

			/*2. Oy direction*/
            std::tie(InviscidTerm[0][1], InviscidTerm[1][1], InviscidTerm[2][1], InviscidTerm[3][1]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, umVal, vVal, vmVal, totalE, pVal, 2);

			return InviscidTerm;
		}

		/*Function calculates Viscous terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, int na, int nb)
		{
            int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

            double uVal(0.0), vVal(0.0), a(0.0), b(0.0);
			/*ViscousTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> ViscousTerm(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> StressHeatFlux(2, std::vector<double>(3, 0.0));
			std::vector<double> vectorU(4, 0.0);
			std::vector<double> vectordUx(4, 0.0);
			std::vector<double> vectordUy(4, 0.0);
			std::tie(a, b) = auxUlti::getGaussCoor(na, nb);

			/*calculate conservative and derivative variables*/
            vectorU[0] = volumeFields::rhoVolGauss[element][nanb];
            vectorU[1] = volumeFields::rhouVolGauss[element][nanb];
            vectorU[2] = volumeFields::rhovVolGauss[element][nanb];
            vectorU[3] = volumeFields::rhoEVolGauss[element][nanb];

            uVal=vectorU[1]/vectorU[0];
            vVal=vectorU[2]/vectorU[0];

            /*
            for (int i = 0; i < 4; i++)
			{
                vectordUx[i] = math::pointAuxValue(element, a, b, i + 1, 1)*muVal;
                vectordUy[i] = math::pointAuxValue(element, a, b, i + 1, 2)*muVal;
            }
            */
            vectordUx=math::pointSVars(0,element,a,b,1,1);
            vectordUy=math::pointSVars(0,element,a,b,2,1);

			/*calculate stresses and heat fluxes*/
            StressHeatFlux = math::viscousTerms::calcStressTensorAndHeatFlux(vectorU, vectordUx, vectordUy);
            std::tie(ViscousTerm[0][0], ViscousTerm[1][0], ViscousTerm[2][0], ViscousTerm[3][0]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 1);
            std::tie(ViscousTerm[0][1], ViscousTerm[1][1], ViscousTerm[2][1], ViscousTerm[3][1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 2);

            //Correct Diffusive Term of NSF Eqns
            //Theo mo hinh Durst
            if (extNSF_Durst::enable)
            {
                //SD=self diffusion
                std::vector<std::vector<double>> diffTerms(2, std::vector<double>(4, 0.0));

                extNSF_Durst::correctViscousTerms(diffTerms,vectorU,vectordUx,vectordUy);
                //Add to classical diffusive terms
                for (int i=0; i<4; i++)
                {
                    ViscousTerm[i][0] = ViscousTerm[i][0] + diffTerms[0][i];
                    ViscousTerm[i][1] = ViscousTerm[i][1] + diffTerms[1][i];
                }
            }

            return ViscousTerm;
		}

        std::tuple<double, double> getInterfacesFluxes(int edge, int element, int nG, int mod, int valType)
        {
            /*
            -mod: 1 for inviscid, 2 for viscous
            -direction: 1 for Ox, 2 for Oy
            */
            bool isMaster(auxUlti::checkMaster(element, edge));
            int locationPlus(-1), locationMinus(-1);
            double valPlus(0.0), valMinus(0.0);
            if (isMaster)
            {
                locationPlus = nG;
                locationMinus = nG + mathVar::nGauss1D + 1;
            }
            else
            {
                locationPlus = nG + mathVar::nGauss1D + 1;
                locationMinus = nG;
            }

            switch (mod)
            {
            case 1: //inviscid
            {
                switch (valType)
                {
                case 1:
                {
                    valPlus = surfaceFields::invis_rho[edge][locationPlus];
                    valMinus = surfaceFields::invis_rho[edge][locationMinus];
                }
                break;
                case 2:
                {
                    valPlus = surfaceFields::invis_rhou[edge][locationPlus];
                    valMinus = surfaceFields::invis_rhou[edge][locationMinus];
                }
                break;
                case 3:
                {
                    valPlus = surfaceFields::invis_rhov[edge][locationPlus];
                    valMinus = surfaceFields::invis_rhov[edge][locationMinus];
                }
                break;
                case 4:
                {
                    valPlus = surfaceFields::invis_rhoE[edge][locationPlus];
                    valMinus = surfaceFields::invis_rhoE[edge][locationMinus];
                }
                break;
                default:
                    break;
                }
            }
                break;
            case 2: //viscous
            {
                switch (valType)
                {
                case 1:
                {
                    valPlus = surfaceFields::Vis_rho[edge][locationPlus];
                    valMinus = surfaceFields::Vis_rho[edge][locationMinus];
                }
                break;
                case 2:
                {
                    valPlus = surfaceFields::Vis_rhou[edge][locationPlus];
                    valMinus = surfaceFields::Vis_rhou[edge][locationMinus];
                }
                break;
                case 3:
                {
                    valPlus = surfaceFields::Vis_rhov[edge][locationPlus];
                    valMinus = surfaceFields::Vis_rhov[edge][locationMinus];
                }
                break;
                case 4:
                {
                    valPlus = surfaceFields::Vis_rhoE[edge][locationPlus];
                    valMinus = surfaceFields::Vis_rhoE[edge][locationMinus];
                }
                break;
                default:
                    break;
                }
            }
                break;
            default:
                break;
            }
            return std::make_tuple(valPlus, valMinus);
        }

		void calcVolumeIntegralTerms(int element, std::vector<double> &VolIntTerm1, std::vector<double> &VolIntTerm2, std::vector<double> &VolIntTerm3, std::vector<double> &VolIntTerm4)
		{
			/*User's guide:
			All input array have form:
			- number of rows: orderElem*/

            //DURST MODEL ------------------------------------------------------------------------------
            extNSF_Durst::needToRemoveDiffTerm=false;
            //DURST MODEL ------------------------------------------------------------------------------

			std::vector<std::vector<double>> ViscousTerms(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> InviscidTerms(4, std::vector<double>(2, 0.0));
			//std::vector<double> VolInt(4, 0.0);

			std::vector<std::vector<double>>
                GsVolX1(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
                GsVolY1(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),

                GsVolX2(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
                GsVolY2(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),

                GsVolX3(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
                GsVolY3(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),

                GsVolX4(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
                GsVolY4(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D+ 1, 0.0));

            for (int na = 0; na <= mathVar::nGauss2D; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
				{
					//std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					/*A INVISCID TERMS*/
					InviscidTerms = NSFEq::calcGaussInviscidTerm(element, na, nb);
					/*A1. Inviscid term on Ox direction*/
                    GsVolX1[na][nb] = InviscidTerms[0][0];
                    GsVolX2[na][nb] = InviscidTerms[1][0];
                    GsVolX3[na][nb] = InviscidTerms[2][0];
                    GsVolX4[na][nb] = InviscidTerms[3][0];
					/*A2. Inviscid term on Oy direction*/
                    GsVolY1[na][nb] = InviscidTerms[0][1];
                    GsVolY2[na][nb] = InviscidTerms[1][1];
                    GsVolY3[na][nb] = InviscidTerms[2][1];
                    GsVolY4[na][nb] = InviscidTerms[3][1];

                    if (flowProperties::viscous)
                    {
                        /*B VISCOUS TERMS*/
                        ViscousTerms = NSFEq::calcGaussViscousTerm(element, na, nb);
                        /*B1. Viscous term on Ox direction*/
                        GsVolX1[na][nb] += ViscousTerms[0][0];
                        GsVolX2[na][nb] += ViscousTerms[1][0];
                        GsVolX3[na][nb] += ViscousTerms[2][0];
                        GsVolX4[na][nb] += ViscousTerms[3][0];
                        /*B2. Viscous term on Oy direction*/
                        GsVolY1[na][nb] += ViscousTerms[0][1];
                        GsVolY2[na][nb] += ViscousTerms[1][1];
                        GsVolY3[na][nb] += ViscousTerms[2][1];
                        GsVolY4[na][nb] += ViscousTerms[3][1];
                    }
				}
			}

            for (int order = 1; order <= mathVar::orderElem; order++)
			{
				/*CALCULATE INTEGRALS*/
                VolIntTerm1[order] = process::volumeInte(element, GsVolX1, order, 1);
                VolIntTerm1[order] += process::volumeInte(element, GsVolY1, order, 2);

                VolIntTerm2[order] = process::volumeInte(element, GsVolX2, order, 1);
                VolIntTerm2[order] += process::volumeInte(element, GsVolY2, order, 2);

                VolIntTerm3[order] = process::volumeInte(element, GsVolX3, order, 1);
                VolIntTerm3[order] += process::volumeInte(element, GsVolY3, order, 2);

                VolIntTerm4[order] = process::volumeInte(element, GsVolX4, order, 1);
                VolIntTerm4[order] += process::volumeInte(element, GsVolY4, order, 2);
			}
            VolIntTerm1[0] = 0;
            VolIntTerm2[0] = 0;
            VolIntTerm3[0] = 0;
            VolIntTerm4[0] = 0;

            //DURST MODEL ------------------------------------------------------------------------------
            //Correct Volume integral if Durst model is enabled
            if (extNSF_Durst::enable)
            {
                extNSF_Durst::correctEnergyEqnVolIntTerm(element,VolIntTerm4);
            }
            //DURST MODEL ------------------------------------------------------------------------------
		}

		void calcSurfaceIntegralTerms(int element, std::vector<double> &SurfIntTerm1, std::vector<double> &SurfIntTerm2, std::vector<double> &SurfIntTerm3, std::vector<double> &SurfIntTerm4)
		{
			/*User's guide:
            All input arrays have to have the form:
			- number of rows: orderElem*/

			int elemType(auxUlti::checkType(element)), edgeName(0);
			int faceBcType(0);

            std::vector<double> Flux1Temp(mathVar::nGauss1D + 1, 0.0),
                Flux2Temp(mathVar::nGauss1D + 1, 0.0),
                Flux3Temp(mathVar::nGauss1D + 1, 0.0),
                Flux4Temp(mathVar::nGauss1D + 1, 0.0);

            std::vector<double> Fluxes(4, 0.0);

			std::vector<std::vector<double>>
                Flux1(mathVar::nGauss1D + 1, std::vector<double>(4, 0.0)),
                Flux2(mathVar::nGauss1D + 1, std::vector<double>(4, 0.0)),
                Flux3(mathVar::nGauss1D + 1, std::vector<double>(4, 0.0)),
                Flux4(mathVar::nGauss1D + 1, std::vector<double>(4, 0.0));

			//std::vector<double> SurInt(4, 0.0);

            for (int nface = 0; nface < elemType; nface++)
			{
                edgeName = meshVar::inedel[element][nface];
				faceBcType = auxUlti::getBCType(edgeName);

                if (faceBcType == 0)  //internal edge
                {
                    for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
					{
						Fluxes = process::NSFEq::getGaussVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss);
                        Flux1[nGauss][nface] = Fluxes[0];
                        Flux2[nGauss][nface] = Fluxes[1];
                        Flux3[nGauss][nface] = Fluxes[2];
                        Flux4[nGauss][nface] = Fluxes[3];
					}
				}
				else  //boundary edge
				{
                    //Tinh he so C neu flux type la LxF
                    if (DGSchemes::fluxControl::LxF)
                    {
                        /* Ham nay lay surfaceFields U o phia - (phia ghost cell) de tinh he so C.
                         * Sau khi chay buoc giai phuong trinh phu, cac bien trong surfaceFields phia ghost cell tai cac bien
                         * da duoc cap nhat theo U tai time ngay lien truoc (n)
                        */
                        BCSupportFncs::NSFEqBCs::findMaxLxFConstantOnEdge(element,edgeName);
                    }

                    for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
					{
						Fluxes = NSFEqBCsImplement(element, edgeName, nGauss);
                        Flux1[nGauss][nface] = Fluxes[0];
                        Flux2[nGauss][nface] = Fluxes[1];
                        Flux3[nGauss][nface] = Fluxes[2];
                        Flux4[nGauss][nface] = Fluxes[3];
					}
				}
			}

            for (int order = 0; order <= mathVar::orderElem; order++)
			{
                for (int nface = 0; nface < elemType; nface++)
				{
                    edgeName = meshVar::inedel[element][nface];  //A BIG BUG!!!!!!
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
					{
                        Flux1Temp[nG] = Flux1[nG][nface];
                        Flux2Temp[nG] = Flux2[nG][nface];
                        Flux3Temp[nG] = Flux3[nG][nface];
                        Flux4Temp[nG] = Flux4[nG][nface];
					}
                    SurfIntTerm1[order] += process::surfaceInte(element, edgeName, Flux1Temp, order);
                    SurfIntTerm2[order] += process::surfaceInte(element, edgeName, Flux2Temp, order);
                    SurfIntTerm3[order] += process::surfaceInte(element, edgeName, Flux3Temp, order);
                    SurfIntTerm4[order] += process::surfaceInte(element, edgeName, Flux4Temp, order);
				}
			}
			//return SurInt;
		}

        std::vector<double> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, int nGauss)
		{
            std::vector<double> Fluxes(4, 0.0);
			double
                invisTerm1P(0.0),
                invisTerm2P(0.0),
                invisTerm3P(0.0),
                invisTerm4P(0.0),
                visTerm1P(0.0),
                visTerm2P(0.0),
                visTerm3P(0.0),
                visTerm4P(0.0);

			/*INVISCID TERM*/
			//Get value
            std::tie(invisTerm1P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 1, 1);
            std::tie(invisTerm2P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 1, 2);
            std::tie(invisTerm3P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 1, 3);
            std::tie(invisTerm4P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 1, 4);
			
			/*VISCOUS TERM*/
			//Get value
            if (flowProperties::viscous)
            {
                std::tie(visTerm1P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 2, 1);
                std::tie(visTerm2P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 2, 2);
                std::tie(visTerm3P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 2, 3);
                std::tie(visTerm4P, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 2, 4);
            }

			/*Calculate fluxes*/
            Fluxes[0] = invisTerm1P+visTerm1P;
            Fluxes[1] = invisTerm2P+visTerm2P;
            Fluxes[2] = invisTerm3P+visTerm3P;
            Fluxes[3] = invisTerm4P+visTerm4P;

			return Fluxes;
		}

		std::vector<double> solveTimeMarching(int element, std::vector<double> &ddtArr, std::vector<double> &UnArr, int RKOrder, int varType)
		{
			std::vector<double> OutArr(mathVar::orderElem + 1, 0.0);
			
			switch (RKOrder)
			{
			case 1:
			{
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					OutArr[order] = dt * ddtArr[order] + UnArr[order];
				}
			}
			break;
			case 2:
			{
				switch (varType)
				{
				case 1://rho
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rho0[element][order];
					}
				}
				break;
				case 2://rhou
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rhou0[element][order];
					}
				}
				break;
				case 3://rhov
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rhov0[element][order];
					}
				}
				break;
				case 4://rhoE
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rhoE0[element][order];
					}
				}
				break;
				default:
					break;
				}
			}
			break;
			case 3:
			{
				switch (varType)
				{
				case 1://rho
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rho0[element][order];
					}
				}
				break;
				case 2://rhou
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rhou0[element][order];
					}
				}
				break;
				case 3://rhov
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rhov0[element][order];
					}
				}
				break;
				case 4://rhoE
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rhoE0[element][order];
					}
				}
				break;
				default:
					break;
				}
			}
			break;
			default:
				break;
			}
			
			return OutArr;
		}
	}

	namespace timeDiscretization
	{
        void solveTimeEq()
        {
            if (systemVar::ddtScheme==1)
            {
                process::timeDiscretization::Euler();
            }
            else if (systemVar::ddtScheme==2)
            {
                //not support
            }
            else if (systemVar::ddtScheme==3)
            {
                process::timeDiscretization::TVDRK3();
            }
        }

		void calcGlobalTimeStep()
		{
            if (systemVar::firstIter)
            {
                dt = 0;
            }
            else
            {
                std::vector<double>timeStepArr(meshVar::nelem2D, 1.0);
                for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
                {
                    timeStepArr[nelement] = process::timeDiscretization::localTimeStep(nelement);

                }
                runTime += dt;
                dt = *std::min_element(timeStepArr.begin(), timeStepArr.end());  //find min value of vector
            }
		}

		double localTimeStep(int element)
		{
			double deltaT(0.0), uVal(0.0), vVal(0.0), velocity(0.0), TVal(0.0), aSound(0.0), muVal(0.0), size(meshVar::cellSize[element]);
			uVal = rhou[element][0] / rho[element][0];
			vVal = rhov[element][0] / rho[element][0];

            /* Dung T convective*/
            TVal = math::CalcTFromConsvVar(rho[element][0], rhou[element][0], rhov[element][0], rhoE[element][0]);
            muVal = math::CalcVisCoef(TVal);

            velocity = sqrt(pow(uVal, 2) + pow(vVal, 2));
			aSound = math::CalcSpeedOfSound(TVal);
			deltaT = fabs((1.0 / pow(mathVar::orderElem + 1, 2))*(size*systemVar::CFL) / (aSound + velocity + muVal/size));
			return deltaT;
		}

		double localErrorEstimate(int element, std::vector<double> &ddtArr)
		{
			double xC(-1.0 / 3.0), yC(1.0 / 3.0), errorVal(0.0);

			if (auxUlti::checkType(element) == 4)
			{
				xC = 0.0;
				yC = 0.0;
			}

			//Compute basis function
			math::basisFc(xC, yC, auxUlti::checkType(element));
            for (int order = 0; order <= mathVar::orderElem; order++)
			{
				errorVal += ddtArr[order] * mathVar::B[order];
			}

			return abs(errorVal);
		}

		void globalErrorEstimate()
		{
			double xC(0.0), yC(0.0), aC(0.0), bC(0.0);
			std::vector<double> rhoRes(meshVar::nelem2D, 0.0),
				rhouRes(meshVar::nelem2D, 0.0),
				rhovRes(meshVar::nelem2D, 0.0),
				rhoERes(meshVar::nelem2D, 0.0);
			double rhoResGlobal(1.0), rhouResGlobal(1.0), rhovResGlobal(1.0), rhoEResGlobal(1.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				std::tie(xC, yC) = auxUlti::getCellCentroid(nelement);
				std::tie(aC, bC) = math::inverseMapping(nelement, xC, yC);
				rhoRes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 1);
				rhouRes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 2);
				rhovRes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 3);
				rhoERes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 4);
			}
			rhoResGlobal = *std::max_element(rhoRes.begin(), rhoRes.end());
			rhouResGlobal = *std::max_element(rhouRes.begin(), rhouRes.end());
			rhovResGlobal = *std::max_element(rhovRes.begin(), rhovRes.end());
			rhoEResGlobal = *std::max_element(rhoERes.begin(), rhoERes.end());

            if (systemVar::currentProc==0)
            {
                std::vector<double> vectorRhoRes(systemVar::totalProc,0.0),
                        vectorRhouRes(systemVar::totalProc,0.0),
                        vectorRhovRes(systemVar::totalProc,0.0),
                        vectorRhoERes(systemVar::totalProc,0.0);

                vectorRhoRes[0]=rhoResGlobal;
                vectorRhouRes[1]=rhouResGlobal;
                vectorRhovRes[2]=rhovResGlobal;
                vectorRhoERes[3]=rhoEResGlobal;

                for (int iproc=1; iproc<systemVar::totalProc; iproc++)
                {
                    MPI_Recv(&vectorRhoRes[iproc], 1, MPI_DOUBLE, iproc, iproc*10+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&vectorRhouRes[iproc], 1, MPI_DOUBLE, iproc, iproc*10+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&vectorRhovRes[iproc], 1, MPI_DOUBLE, iproc, iproc*10+3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&vectorRhoERes[iproc], 1, MPI_DOUBLE, iproc, iproc*10+4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                IO::residualOutput(*std::max_element(vectorRhoRes.begin(), vectorRhoRes.end()),
                                   *std::max_element(vectorRhouRes.begin(), vectorRhouRes.end()),
                                   *std::max_element(vectorRhovRes.begin(), vectorRhovRes.end()),
                                   *std::max_element(vectorRhoERes.begin(), vectorRhoERes.end()));
            }
            else
            {
                MPI_Send(&rhoResGlobal, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+1, MPI_COMM_WORLD);
                MPI_Send(&rhouResGlobal, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+2, MPI_COMM_WORLD);
                MPI_Send(&rhovResGlobal, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+3, MPI_COMM_WORLD);
                MPI_Send(&rhoEResGlobal, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+4, MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            if (systemVar::iterCount==1)
            {
                if (systemVar::currentProc==0)
                {
                    //Send systemVar::ResNorm den cac proc khac
                    for (int iproc=1; iproc<systemVar::totalProc; iproc++)
                    {
                        MPI_Send(&systemVar::rhoResNorm, 1, MPI_DOUBLE, iproc, iproc*10+1, MPI_COMM_WORLD);
                        MPI_Send(&systemVar::rhouResNorm, 1, MPI_DOUBLE, iproc, iproc*10+2, MPI_COMM_WORLD);
                        MPI_Send(&systemVar::rhovResNorm, 1, MPI_DOUBLE, iproc, iproc*10+3, MPI_COMM_WORLD);
                        MPI_Send(&systemVar::rhoEResNorm, 1, MPI_DOUBLE, iproc, iproc*10+4, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Recv(&systemVar::rhoResNorm, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&systemVar::rhouResNorm, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&systemVar::rhovResNorm, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&systemVar::rhoEResNorm, 1, MPI_DOUBLE, 0, systemVar::currentProc*10+4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
		}

        /**
         * @brief Function solves 1 TVD_RK iteration
         *
         * Workflow explanation
         *
         * Assume that U is at time n:
         * - 1: Run Function Positivity preserving limiter for Rho (to calculate theta1)
         * to keep Rho always positive while solving \f$S_m = \nabla (\rho)\f$ in case of
         * mass diffusion = ON. If mass diffusion = OFF, Limiter
         * for all remaining conservative variable will be call right after limiter function
         * for Rho to limit vector U.
         *
         * - 2: If mass diffusion = ON, calculate value Rho at all volume and surface Gauss points
         * and then solve \f$S_m = \nabla (\rho)\f$ equation. Results are saved into BR1Vars::massDiffusion::rhoX
         * and BR1Vars::massDiffusion::rhoY.
         *
         * - 3: Apply Limiter for vector U.
         *
         * - 4: Solve T from U.
         *
         * - 5: Solve equation \f$S = \mu \nabla (\rho)\f$.
         *
         * - 6: Calculate fluxes (inviscid va viscous) at all faces.
         *
         * - 7: Solve NSF
         *
         * @param RKOrder: Runge-Kutta order
         */
		void TVDRK_1step(int RKOrder)
		{
            limiter::limiter_1InnerStep();

            //COMPUTE GAUSS VALUES AND INTERFACES VALUES
			process::calcVolumeGaussValues();
            process::calcValuesAtInterface();

            //CALCULATE T
            process::calcTGauss();

            //SOLVE AUXILARY EQUATION
            if (flowProperties::viscous)
            {
                process::auxEq::solveAuxEquation();
            }
            else
            {
                //Chay ham nay de update gia tri U tai surfaceBCFields
                process::auxEq::updateSurfaceFieldsAtBC();
            }

            //Show warning
            message::showWarning();

			//SOLVE NSF EQUATION
            process::NSFEq::calcSurfaceFlux();
            process::NSFEq::solveNSFEquation_TVDRK(RKOrder);
		}

        /**
         * @brief Function solves 1 time steps of DG workflow TVDRK (otal Variation Diminishing Runge-Kutta) time discretization
         */
		void TVDRK3()
        {
            rho0 = rho;
            rhou0 = rhou;
            rhov0 = rhov;
            rhoE0 = rhoE;

            for (int iRKOrder = 1; iRKOrder <= 3; iRKOrder++)
			{
				process::timeDiscretization::TVDRK_1step(iRKOrder);
                rho = rhoN;
                rhou = rhouN;
                rhov = rhovN;
                rhoE = rhoEN;
			}

            //UPDATE TIME VARYING BC
            timeVaryingBCs::updateBCs();
            //UPDATE NONEQUILIBRIUM BC
            nonEquilibriumBCs::updateBCs();

            limiter::limiter_1OutterStep();
		}

        void Euler()
        {
            rho0 = rho;
            rhou0 = rhou;
            rhov0 = rhov;
            rhoE0 = rhoE;
            //---------CALCULATION OF EACH EULER ITERATION---------//
            limiter::limiter_1InnerStep();

            //COMPUTE GAUSS VALUES AND INTERFACES VALUES
            process::calcVolumeGaussValues();
            process::calcValuesAtInterface();

            //CALCULATE T
            process::calcTGauss();

            //SOLVE AUXILARY EQUATION
            if (flowProperties::viscous)
            {
                process::auxEq::solveAuxEquation();
            }
            else
            {
                //Chay ham nay de update gia tri U tai surfaceBCFields
                process::auxEq::updateSurfaceFieldsAtBC();
            }

            //Show warning
            message::showWarning();

            //SOLVE NSF EQUATION
            process::NSFEq::calcSurfaceFlux();
            process::NSFEq::solveNSFEquation_Euler();
            //---------CALCULATION OF EACH EULER ITERATION---------//

            rho = rhoN;
            rhou = rhouN;
            rhov = rhovN;
            rhoE = rhoEN;

            //UPDATE TIME VARYING BC
            timeVaryingBCs::updateBCs();
            //UPDATE NONEQUILIBRIUM BC
            nonEquilibriumBCs::updateBCs();

            limiter::limiter_1OutterStep();
        }

        namespace parallel {
        void minTimeStep()
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (systemVar::currentProc==0)
            {
                std::vector<double>vectorTimeStep(systemVar::totalProc,1.0);
                vectorTimeStep[0]=dt;
                for (int inode=1;inode<systemVar::totalProc;inode++)
                {
                    MPI_Recv(&vectorTimeStep[inode], 1, MPI_DOUBLE, inode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                dt=*std::min_element(vectorTimeStep.begin(), vectorTimeStep.end());
            }
            else
            {
                MPI_Send(&dt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }

            //Sau khi tim duoc dt min, send gia tri dt lai ve cac process
            if (systemVar::currentProc==0)
            {
                for (int inode=1;inode<systemVar::totalProc;inode++)
                {
                    MPI_Send(&dt, 1, MPI_DOUBLE, inode, 0, MPI_COMM_WORLD);
                }
            }
            else
            {
                MPI_Recv(&dt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        }
	}

	std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType)
	{
		//valType = 1(rho*mu), 2(rhou*mu), 3(rhov*mu), 4(rhoE*mu)
		bool isMaster(auxUlti::checkMaster(element, edge));
		int locationPlus(-1), locationMinus(-1);
		double valPlus(0.0), valMinus(0.0);
		if (isMaster)
		{
			locationPlus = nG;
            locationMinus = nG + mathVar::nGauss1D + 1;
		}
		else
		{
            locationPlus = nG + mathVar::nGauss1D + 1;
			locationMinus = nG;
		}

		switch (valType)
		{
		case 1:
		{
            valPlus = surfaceFields::rho[edge][locationPlus];
            valMinus = surfaceFields::rho[edge][locationMinus];
		}
		break;
		case 2:
		{
            valPlus = surfaceFields::rhou[edge][locationPlus];
            valMinus = surfaceFields::rhou[edge][locationMinus];
		}
		break;
		case 3:
		{
            valPlus = surfaceFields::rhov[edge][locationPlus];
            valMinus = surfaceFields::rhov[edge][locationMinus];
		}
		break;
		case 4:
		{
            valPlus = surfaceFields::rhoE[edge][locationPlus];
            valMinus = surfaceFields::rhoE[edge][locationMinus];
		}
		break;
		default:
			break;
		}
		return std::make_tuple(valPlus, valMinus);
	}

	double volumeInte(int elem, std::vector< std::vector<double> > &Ui, int order, int direction)
	{
		double Int(0.0);
        std::vector<std::vector<double>> A(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));
		double dBi(0.0);

        for (int na = 0; na <= mathVar::nGauss2D; na++)
		{
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
			{
				dBi = math::Calc_dBxdBy(elem, order, na, nb, direction);
                A[na][nb] = dBi * Ui[na][nb];
			}
		}
		Int = math::volumeInte(A, elem);

		return Int;
	}

	double surfaceInte(int elem, int edge, std::vector<double> &FluxVector, int order)
	{
		double Bi(0.0), a(0.0), b(0.0), inte(0.0);
        std::vector<double> F(mathVar::nGauss1D + 1, 0.0);

		//inte = 0.0;  //value of integral is reset for each order
        for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
		{
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, elem, nG);
			math::basisFc(a, b, auxUlti::checkType(elem));  //mathVar::B is changed after this command line is excuted
			Bi = mathVar::B[order];
			F[nG] = Bi * FluxVector[nG];
		}
        inte = math::surfaceInte(F, edge);

		return inte;
	}

	std::vector<std::vector<double>> calculateStiffMatrix(int element)
	{
		std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));

        for (int order1 = 0; order1 <= mathVar::orderElem; order1++)
		{
            for (int order2 = 0; order2 <= mathVar::orderElem; order2++)
			{
				StiffMatrix[order1][order2] = process::calculateStiffMatrixElement(element, order1, order2);
			}
		}
		return StiffMatrix;
	}

	double calculateStiffMatrixElement(int element, int order1, int order2)
	{
		double B1(0.0), B2(0.0), Inte(0.0);
        std::vector<std::vector<double>> FMatrix(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));
		int elemType(auxUlti::checkType(element));

		switch (elemType)
		{
		case 3:
		{
            for (int na = 0; na <= mathVar::nGauss2D; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
				{
                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                    B1 = mathVar::BPts_Tri[order1][nanb];
                    B2 = mathVar::BPts_Tri[order2][nanb];
					FMatrix[na][nb] = B1 * B2;
				}
			}
		}
		break;
		case 4:
		{
            for (int na = 0; na <= mathVar::nGauss2D; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
				{
                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                    B1 = mathVar::BPts_Quad[order1][nanb];
                    B2 = mathVar::BPts_Quad[order2][nanb];
					FMatrix[na][nb] = B1 * B2;
				}
			}
		}
		break;
		default:
			break;
		}
		Inte = math::volumeInte(FMatrix, element);
		return Inte;
	}

	bool checkRunningCond()
	{
		bool run(true);
		if (runTime >= systemVar::Ttime)
		{
			run = false;  //stop running if runtime bigger than Ttime (total time)
		}
		return run;
	}
}
