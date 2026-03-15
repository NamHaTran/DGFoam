#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include "DGMessagesLib.h"
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "DGIOLib.h"

//Extended Navier-Stokes-Fourier model
#include "./extNSFEqns/FranzDurst/DurstModel.h"

namespace debugTool
{
	void checkElemSurPt(int ipoin)
	{
		std::cout << "In SALOME, point ";
		std::cout << ipoin << " is surrounded by elements \n";
		ipoin--;
		for (int ii = meshVar::esup2[ipoin]+1; ii <= meshVar::esup2[ipoin+1]; ii++)
		{
			std::cout << meshVar::esup1[ii]+ meshVar::nelem1D + 1 << std::endl;
		}
	}

	void checkPtsSurPt(int ipoin)
	{
		std::cout << "In SALOME, point ";
		std::cout << ipoin << " is surrounded by points \n";
		ipoin--;
		for (int ii = meshVar::psup2[ipoin] + 1; ii <= meshVar::psup2[ipoin + 1]; ii++)
		{
			std::cout << meshVar::psup1[ii] + 1 << std::endl;
		}
	}

	void checkElemsSurElem(int ielem)
	{
		std::cout << "In SALOME, element ";
		std::cout << ielem << " is surrounded by elements \n";
		ielem = ielem - 1 - meshVar::nelem1D;
		for (int i = 0; i < 4; i++)
		{
            if (meshVar::esuel[ielem][i]>=0)
			{
                std::cout << meshVar::esuel[ielem][i] + 1 + meshVar::nelem1D << std::endl;
			}
		}
	}

	void checkElemInfor(int elem)
	{
		std::cout << "In SALOME, element ";
		std::cout << elem << " has the following edges \n";
		elem = elem - 1 - meshVar::nelem1D;
		int edgeName(0), point1(0), point2(0), grp(0), bcType(0), edgeOrder(0), elemType(0);
		double nx(0.0), ny(0.0);
		bool master(true);
		elemType = auxUlti::checkType(elem);
		
		if (elem>0)
		{
			for (int iedge = 0; iedge < elemType; iedge++)
			{
                edgeName = meshVar::inedel[elem][iedge];
                point1 = meshVar::inpoed[edgeName][0] + 1;
                point2 = meshVar::inpoed[edgeName][1] + 1;
				grp = auxUlti::getGrpOfEdge(edgeName);
				bcType = auxUlti::getBCType(edgeName);
				nx = (auxUlti::getNormVectorComp(elem, edgeName, 1));
				ny = (auxUlti::getNormVectorComp(elem, edgeName, 2));
				edgeOrder = auxUlti::findEdgeOrder(elem, edgeName);
				master = auxUlti::checkMaster(elem, edgeName);
				std::cout << "+ Edge " << edgeName << " created by two points: " << point1 << " and " << point2 << std::endl
					<< "	Coordinates of normal vector are: " << nx << ", " << ny << std::endl
					<< "	Edge belongs to group " << grp << " and has bc type is " << bcType << std::endl
					<< "	Order of edge is " << edgeOrder << std::endl
                    << "	Edge is shared by two elements that are " << meshVar::ineled[edgeName][0] + 1 + meshVar::nelem1D << " and " << meshVar::ineled[edgeName][1] + 1 + meshVar::nelem1D << std::endl;
				if (master)
				{
					std::cout << "	Considering element is a master of edge\n";
				}
				else
				{
					std::cout << "	Considering element is not a master of edge\n";
				}
			}
		}
		else
		{
			std::cout << "	Considering element is not a polygonal element!\n";
		}

	}

	void checkPointValue(int element)
	{
		double rhoVal(math::pointValue(element, 0.0, 0.0, 1, 1)),
			uVal(math::pointValue(element, 0.0, 0.0, 2, 1)),
			vVal(math::pointValue(element, 0.0, 0.0, 3, 1)),
			//eVal(math::pointValue(element, 0.0, 0.0, 4, 1)),
			pVal(math::pointValue(element, 0.0, 0.0, 5, 1)),
			TVal(math::pointValue(element, 0.0, 0.0, 6, 1));
			//muVal(math::pointValue(element, 0.0, 0.0, 7, 1));
		std::cout << "- Element: " << element + meshVar::nelem1D + 1 << std::endl
			<< "- rho: " << rhoVal << std::endl
			<< "- u: " << uVal << std::endl
			<< "- v: " << vVal << std::endl
			//<< "- e: " << eVal << std::endl
			<< "- p: " << pVal << std::endl
			<< "- T: " << TVal << std::endl << std::endl;
			//<< "- mu: " << muVal << std::endl << std::endl;
	}
}

void reconstructLatestTime()
{
    /* Chu y!!!:
     * - Chi reconstruct duoc cac fields rho, rhou, rhov, rhoE.
     * - Chua reconstruct duoc surfaceBCFields
     * - Khi chay lenh reconstructcase, phai setup cac dieu kien trong file DGOption
     * ve time dang muon reconstruct, bao gom:
     *  + totalProcess
     *  + orderOfAccuracy
     * - Set time trong file time ve time muon reconstruct
     * - Chua su dung duoc cho truong hop mass diffusion -> khi mass diffusion on, can save
     * data div(U) xuong file.
    */

    //Resize array
    rho = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    rhou = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    rhov = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    rhoE = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

    SurfaceBCFields::uBc = auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss1D, 0.0);
    SurfaceBCFields::vBc = auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss1D, 0.0);
    SurfaceBCFields::TBc = auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss1D, iniValues::TIni);
    SurfaceBCFields::pBc = auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss1D, iniValues::pIni);

    auxUlti::initialize1DArray(theta1Arr, meshVar::nelem2D, 1.0);
    theta2Arr=new double[meshVar::nelem2D];
    auxUlti::initialize1DArray(theta2Arr, meshVar::nelem2D, 1.0);
    BR1Vars::rhoX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    BR1Vars::rhouX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    BR1Vars::rhovX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    BR1Vars::rhoEX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

    BR1Vars::rhoY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    BR1Vars::rhouY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    BR1Vars::rhovY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    BR1Vars::rhoEY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
    //-------------------------------------------

    int *rankOf2DElem;
    int *Elem2DlocalIdWithRank;
    int startId, endId, a, temp;

    std::vector<std::vector<double>> rhoTemp(meshVar::nelem2D,std::vector<double>(mathVar::orderElem+1,0.0)),
            rhouTemp(meshVar::nelem2D,std::vector<double>(mathVar::orderElem+1,0.0)),
            rhovTemp(meshVar::nelem2D,std::vector<double>(mathVar::orderElem+1,0.0)),
            rhoETemp(meshVar::nelem2D,std::vector<double>(mathVar::orderElem+1,0.0));

    std::string rankOf2DElemLoc = systemVar::pwd + "/Constant/Mesh/rankOf2DElem.txt",
            Elem2DlocalIdWithRankLoc = systemVar::pwd + "/Constant/Mesh/Elem2DlocalIdWithRank.txt";
    std::tie(std::ignore,rankOf2DElem,temp)=IO::read1DIntArray(rankOf2DElemLoc,"rankOf2DElem.txt",true);
    std::tie(std::ignore,Elem2DlocalIdWithRank,temp)=IO::read1DIntArray(Elem2DlocalIdWithRankLoc,"Elem2DlocalIdWithRank.txt",true);

    std::vector<int> endIdList(systemVar::totalProc,0);
    std::string rhoLoc, rhouLoc, rhovLoc, rhoELoc;

    //Doc ket qua tu tat ca cac processor va luu vao cac array *temp
    for (int iproc=0; iproc<systemVar::totalProc; iproc++)
    {
        //Doc cac file ket qua tu cac processor
        double **rhoProc = new double *[1];
        double **rhouProc = new double *[1];
        double **rhovProc = new double *[1];
        double **rhoEProc = new double *[1];

        rhoLoc = systemVar::pwd + "/Processor" + std::to_string(iproc) + "/" + std::to_string(systemVar::iterCount) + "/rho.txt";
        rhouLoc = systemVar::pwd + "/Processor" + std::to_string(iproc) + "/" + std::to_string(systemVar::iterCount) + "/rhou.txt";
        rhovLoc = systemVar::pwd + "/Processor" + std::to_string(iproc) + "/" + std::to_string(systemVar::iterCount) + "/rhov.txt";
        rhoELoc = systemVar::pwd + "/Processor" + std::to_string(iproc) + "/" + std::to_string(systemVar::iterCount) + "/rhoE.txt";

        std::tie(std::ignore,rhoProc,a) = IO::read2DArray(mathVar::orderElem+1,rhoLoc,"rho.txt",true);
        if (iproc==0)
        {
            endIdList[iproc]=a;
        }
        else
        {
            endIdList[iproc]=endIdList[iproc-1]+a;
        }
        std::tie(std::ignore,rhouProc,a)=IO::read2DArray(mathVar::orderElem+1,rhouLoc,"rhou.txt",true);
        std::tie(std::ignore,rhovProc,a)=IO::read2DArray(mathVar::orderElem+1,rhovLoc,"rhov.txt",true);
        std::tie(std::ignore,rhoEProc,a)=IO::read2DArray(mathVar::orderElem+1,rhoELoc,"rhoE.txt",true);

        if (iproc==0)
        {
            startId=0;
            endId=endIdList[0]-1;
        }
        else
        {
            startId=endIdList[iproc-1];
            endId=endIdList[iproc]-1;
        }
        for (int i=startId; i<=endId; i++)
        {
            for (int iorder=0; iorder<mathVar::orderElem+1; iorder++)
            {
                rhoTemp[i][iorder]=rhoProc[i-startId][iorder];
                rhouTemp[i][iorder]=rhouProc[i-startId][iorder];
                rhovTemp[i][iorder]=rhovProc[i-startId][iorder];
                rhoETemp[i][iorder]=rhoEProc[i-startId][iorder];
            }
        }
    }

    //Reconstruct cac field
    std::cout<<"Reconstructing fields...\n" << std::endl;
    int localId, IdOnTempArr, proc;
    for (int ielem=0; ielem<meshVar::nelem2D; ielem++)
    {
        proc=rankOf2DElem[ielem];
        if (proc==0)
        {
            startId=0;
        }
        else
        {
            startId=endIdList[proc-1];
        }
        localId=Elem2DlocalIdWithRank[ielem];
        IdOnTempArr=startId+localId;
        for (int iorder=0; iorder<mathVar::orderElem+1; iorder++)
        {
            rho[ielem][iorder]=rhoTemp[IdOnTempArr][iorder];
            rhou[ielem][iorder]=rhouTemp[IdOnTempArr][iorder];
            rhov[ielem][iorder]=rhovTemp[IdOnTempArr][iorder];
            rhoE[ielem][iorder]=rhoETemp[IdOnTempArr][iorder];
        }
    }

    systemVar::parallelMode=false;
    std::cout << "Saving case...\n" << std::endl;
    IO::saveCase_reconstruct();
    std::cout << "Exporting data to Tecplot...\n" << std::endl;
    DG2Tecplot::exportCellCenteredData(systemVar::iterCount);
}

namespace DG2Tecplot
{
	std::vector<double> calcNodeValues(int valType)
	{
		//Only for primary variables
		std::vector<int>elemSurPt;
		std::vector<double>nodeValues(meshVar::npoin,0.0), U(4, 0.0);
		int numElemSurPt(0), elemId(-1), ptAtBCId(-2);
		double aP(0.0), bP(0.0), localValue(0.0);
		for (int ipoint = 0; ipoint < meshVar::npoin; ipoint++)
		{
			ptAtBCId = meshVar::markPointsAtBC[ipoint] - 1;
			if (ptAtBCId < 0)
			{
				localValue = 0.0;
				elemSurPt = auxUlti::postProcess::getElementsSurroundingPoint(ipoint);
                numElemSurPt = static_cast<int>(elemSurPt.size());
				for (int nelem = 0; nelem < numElemSurPt; nelem++)
				{
					elemId = elemSurPt[nelem];
					std::tie(aP, bP) = auxUlti::postProcess::findPointCoorInStandardSpace(ipoint, elemId);
					localValue += math::pointValue(elemId, aP, bP, valType, 1);
				}
				nodeValues[ipoint] = (localValue / numElemSurPt);
			}
			else
			{
				U = DG2Tecplot::calcNodeValuesAtBC(ptAtBCId);
				if (valType == 1)  //rho
				{
					nodeValues[ipoint] = U[0];
				}
				else if (valType == 2)  //u
				{
					double rhoVal(U[0]), rhouVal(U[1]);
					nodeValues[ipoint] = rhouVal / rhoVal;
				}
				else if (valType == 3)  //v
				{
					double rhoVal(U[0]), rhovVal(U[2]);
					nodeValues[ipoint] = rhovVal / rhoVal;
				}
				else if (valType == 4)  //e
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					nodeValues[ipoint] = material::Cv*math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
				}
				else if (valType == 5)  //p
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
					nodeValues[ipoint] = math::CalcP(TVal, rhoVal);
				}
				else if (valType == 6)  //T
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					nodeValues[ipoint] = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
				}
				else if (valType == 7)  //mu
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
					if ((TVal < 0) && fabs(TVal) < 0.001)
					{
						TVal = fabs(TVal);
					}
					nodeValues[ipoint] = math::CalcVisCoef(TVal);
				}
			}
		}
		return nodeValues;
	}

	std::vector<double> calcNodeValuesAtBC(int ptAtBCId)
	{
		std::vector<double> U(4, 0.0);
		int BCEdgeId(0), BCElem(-1);
		double a(0.0), b(0.0);
		for (int iedge = 0; iedge < 2; iedge++)
		{
			BCEdgeId = SurfaceBCFields::BCPointsInfor[ptAtBCId][iedge];
			BCElem = meshVar::MasterElemOfEdge[BCEdgeId];
			int edgeGrp(auxUlti::getGrpOfEdge(BCEdgeId));
            int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);

			std::tie(a, b) = auxUlti::postProcess::findPointCoorInStandardSpace(BCEdgeId, BCElem);

			if (UType == 1 && TType == 1 && pType == 1)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::patch::inFlow(BCElem, edgeGrp, a, b);
			}
			else if (UType == 4 && TType == 4 && pType == 4)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::patch::outFlow(BCElem, edgeGrp, a, b);
			}
			else if (UType == 2 && TType == 2 && pType == 2)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::wall::noSlipIsoThermal(BCElem, edgeGrp, a, b);
			}
			else if (UType == 2 && TType == 3 && pType == 2)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::wall::noSlipAdiabatic(BCElem, a, b);
			}
			else if (UType == 7 && TType == 7 && pType == 7)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::Symmetry(BCElem, BCEdgeId, a, b);
			}
			else
			{
				std::string errorStr = message::BcCompatibleError(edgeGrp);
				message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
			}
		}
		return U;
	}

	std::vector<double> calcCellCenteredValues(int valType)
	{
		std::vector<double>cellCenteredValues(meshVar::nelem2D, 0.0);
        //double xC(0.0), yC(0.0), aC(0.0), bC(0.0);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
            //std::tie(xC, yC) = auxUlti::getCellCentroid(ielem);
            //std::tie(aC, bC) = math::inverseMapping(ielem, xC, yC);
            cellCenteredValues[ielem] = DG2Tecplot::pointMeanValue(ielem, valType);
		}
		return cellCenteredValues;
	}

    double pointMeanValue(int element, int valType)
    {
        //Compute primary variables from conservative variables
        double out(0.0);
        if (valType == 1)  //rho
        {
            out= rho[element][0];
        }
        else if (valType == 2)  //u
        {
            double rhoVal(rho[element][0]),
                rhouVal(rhou[element][0]);

            out=(rhouVal / rhoVal); //convective velocity

            //DURST MODEL ------------------------------------------------------------------------------
            //Correct Velocity if Durst model is enabled
            if (extNSF_Durst::enable)
            {
                double vD(0.0);
                std::vector<double> U(4, 0.0), dU(4, 0.0);
                U[0]=rho[element][0];
                U[1]=rhou[element][0];
                U[2]=rhov[element][0];
                U[3]=rhoE[element][0];

                dU[0]=BR1Vars::rhoX[element][0];
                dU[1]=BR1Vars::rhouX[element][0];
                dU[2]=BR1Vars::rhovX[element][0];
                dU[3]=BR1Vars::rhoEX[element][0];

                vD=extNSF_Durst::calcDiffVelocity(U,dU);

                //Corect velocity
                out = out + vD;
            }
            //DURST MODEL ------------------------------------------------------------------------------
        }
        else if (valType == 3)  //v
        {
            double rhoVal(rho[element][0]),
                rhovVal(rhov[element][0]);

            out = (rhovVal / rhoVal); //convective velocity

            //DURST MODEL ------------------------------------------------------------------------------
            //Correct Velocity if Durst model is enabled
            if (extNSF_Durst::enable)
            {
                double vD(0.0);
                std::vector<double> U(4, 0.0), dU(4, 0.0);
                U[0]=rho[element][0];
                U[1]=rhou[element][0];
                U[2]=rhov[element][0];
                U[3]=rhoE[element][0];

                dU[0]=BR1Vars::rhoY[element][0];
                dU[1]=BR1Vars::rhouY[element][0];
                dU[2]=BR1Vars::rhovY[element][0];
                dU[3]=BR1Vars::rhoEY[element][0];

                vD=extNSF_Durst::calcDiffVelocity(U,dU);

                //Corect velocity
                out = out + vD;
            }
            //DURST MODEL ------------------------------------------------------------------------------
        }
        else if (valType == 4)  //e
        {
            double rhoVal(rho[element][0]),
                rhouVal(rhou[element][0]),
                rhovVal(rhov[element][0]),
                rhoEVal(rhoE[element][0]);
            out = material::Cv*math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
        }
        else if (valType == 5)  //p
        {
            double rhoVal(rho[element][0]),
                rhouVal(rhou[element][0]),
                rhovVal(rhov[element][0]),
                rhoEVal(rhoE[element][0]);
            out = math::CalcP(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal), rhoVal);
        }
        else if (valType == 6)  //T
        {
            double rhoVal(rho[element][0]),
                rhouVal(rhou[element][0]),
                rhovVal(rhov[element][0]),
                rhoEVal(rhoE[element][0]);
            out = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
        }
        else if (valType == 7)  //mu
        {
            double rhoVal(rho[element][0]),
                rhouVal(rhou[element][0]),
                rhovVal(rhov[element][0]),
                rhoEVal(rhoE[element][0]);
            double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
            out = math::CalcVisCoef(TVal);
            if (out < 0 || out != out)
            {
                std::cout << "unphysical mu at cell " << element + meshVar::nelem1D + 1 << std::endl;
                exit(1);
            }
        }
        else if (valType == 10)  //rho
        {
            double rhoVal(rho[element][0]),
                rhouVal(rhou[element][0]),
                rhovVal(rhov[element][0]),
                rhoEVal(rhoE[element][0]);
            double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));

            out= BR1Vars::rhoX[element][0]/math::CalcVisCoef(TVal);
        }
        else if (valType == 11)  //rho
        {
            double rhoVal(rho[element][0]),
                rhouVal(rhou[element][0]),
                rhovVal(rhov[element][0]),
                rhoEVal(rhoE[element][0]);
            double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));

            out= BR1Vars::rhoY[element][0]/math::CalcVisCoef(TVal);
        }
        return out;
    }

	void exportNodeData(int iter)
	{
		std::vector<double>nodeRho(meshVar::npoin, 0.0), node_u(meshVar::npoin, 0.0), node_v(meshVar::npoin, 0.0), node_p(meshVar::npoin, 0.0), node_T(meshVar::npoin, 0.0), node_uMag(meshVar::npoin, 0.0);
		nodeRho = DG2Tecplot::calcNodeValues(1);
		node_u = DG2Tecplot::calcNodeValues(2);
		node_v = DG2Tecplot::calcNodeValues(3);
		node_p = DG2Tecplot::calcNodeValues(5);
		node_T = DG2Tecplot::calcNodeValues(6);
		for (int i = 0; i < meshVar::npoin; i++)
		{
			node_uMag[i] = sqrt(pow(node_u[i], 2) + pow(node_v[i], 2));
		}

		std::string iter_str = std::to_string(iter);
        std::string fileName(systemVar::caseName + "_node.dat"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/TecplotFile/" + iter_str), code;
        auxUlti::createFolder(Loc, true);
		
        std::string fileLoc(Loc + "/" + fileName);
		std::ofstream fileFlux(fileLoc.c_str());

		code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X", "Y", "RHO", "U" ,"V", "VELOCITY_MAG", "P", "T"
ZONE T="ZONE 1"
ZONETYPE=FEQUADRILATERAL
DATAPACKING=BLOCK)";

		if (fileFlux)
		{
			fileFlux << code << std::endl << "NODES=" << std::to_string(meshVar::npoin) << ", " << "ELEMENTS=" << std::to_string(meshVar::nelem2D) << std::endl;
			//X
			int counter(0);
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][0] << " ";
			}
			fileFlux << std::endl;

			//Y
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][1] << " ";
			}
			fileFlux << std::endl;

			//RHO
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << nodeRho[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_u[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_v[i] << " ";
			}
			fileFlux << std::endl;

			//VELOCITY_MAG
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_uMag[i] << " ";
			}
			fileFlux << std::endl;

			//P
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_p[i] << " ";
			}
			fileFlux << std::endl;

			//T
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_T[i] << " ";
			}
			fileFlux << std::endl;

			//CONNECTIVITY
			for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
			{
				int elemType(auxUlti::checkType(ielem));
				for (int ipoin = 0; ipoin < 3; ipoin++)
				{
					fileFlux << meshVar::Elements2D[ielem][ipoin] + 1 << " ";
				}
				switch (elemType)
				{
				case 3:
				{
					fileFlux << meshVar::Elements2D[ielem][0] + 1 << std::endl;
					break;
				}
				case 4:
				{
					fileFlux << meshVar::Elements2D[ielem][3] + 1 << std::endl;
					break;
				}
				default:
					break;
				}
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(systemVar::caseName + ".dat", fileLoc));
		}
	}

    void exportCellCenteredData(int iter)
	{
        std::vector<double>nodeRho(meshVar::nelem2D, 0.0), node_u(meshVar::nelem2D, 0.0), node_v(meshVar::nelem2D, 0.0), node_p(meshVar::nelem2D, 0.0), node_T(meshVar::nelem2D, 0.0), node_uMag(meshVar::nelem2D, 0.0),
                nodeRhoX(meshVar::nelem2D, 0.0), nodeRhoY(meshVar::nelem2D, 0.0);

		nodeRho = DG2Tecplot::calcCellCenteredValues(1);
		node_u = DG2Tecplot::calcCellCenteredValues(2);
		node_v = DG2Tecplot::calcCellCenteredValues(3);
		node_p = DG2Tecplot::calcCellCenteredValues(5);
		node_T = DG2Tecplot::calcCellCenteredValues(6);

        if (flowProperties::viscous)
        {
            nodeRhoX = DG2Tecplot::calcCellCenteredValues(10);
            nodeRhoY = DG2Tecplot::calcCellCenteredValues(11);
        }

		for (int i = 0; i < meshVar::nelem2D; i++)
		{
			node_uMag[i] = sqrt(pow(node_u[i], 2) + pow(node_v[i], 2));
		}

		std::string iter_str = std::to_string(iter);
        std::string fileName(systemVar::caseName + "cellCentered.dat"), Loc, code;
        Loc = auxUlti::createTimeStepFolder(iter,"tecplot");

        std::string fileLoc(Loc + "/" + fileName);
		std::ofstream fileFlux(fileLoc.c_str());

		code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X"
"Y"
"Rho"
"u"
"v"
"Velocity_mag"
"P"
"T"
"dRhoX"
"dRhoY"
"Theta1"
"Theta2"
ZONE T="ZONE 1"
STRANDID=0
ZONETYPE=FEQuadrilateral
DATAPACKING=BLOCK
VARLOCATION=([3-12]=CELLCENTERED)
DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)
)";

		if (fileFlux)
		{
			fileFlux << code << std::endl << "Nodes=" << std::to_string(meshVar::npoin) << ", " << "Elements=" << std::to_string(meshVar::nelem2D) << std::endl;
			//X
			int counter(0);
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][0] << " ";
			}
			fileFlux << std::endl;

			//Y
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][1] << " ";
			}
			fileFlux << std::endl;

			//RHO
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << nodeRho[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_u[i] << " ";
			}
			fileFlux << std::endl;

			//V
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_v[i] << " ";
			}
			fileFlux << std::endl;

			//VELOCITY_MAG
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_uMag[i] << " ";
			}
			fileFlux << std::endl;

			//P
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_p[i] << " ";
			}
			fileFlux << std::endl;

			//T
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_T[i] << " ";
			}
			fileFlux << std::endl;

			//Theta1
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
                fileFlux << nodeRhoX[i] << " ";
			}
			fileFlux << std::endl;

			//Theta2
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
                fileFlux << nodeRhoY[i] << " ";
			}
			fileFlux << std::endl;

            //Theta1
            counter = 0;
            for (int i = 0; i < meshVar::nelem2D; i++)
            {
                counter++;
                if (counter == 5)
                {
                    fileFlux << std::endl;
                    counter = 0;
                }
                fileFlux << theta1Arr[i] << " ";
            }
            fileFlux << std::endl;

            //Theta2
            counter = 0;
            for (int i = 0; i < meshVar::nelem2D; i++)
            {
                counter++;
                if (counter == 5)
                {
                    fileFlux << std::endl;
                    counter = 0;
                }
                fileFlux << theta2Arr[i] << " ";
            }
            fileFlux << std::endl;

			//CONNECTIVITY
			for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
			{
				int elemType(auxUlti::checkType(ielem));
				for (int ipoin = 0; ipoin < 3; ipoin++)
				{
					fileFlux << meshVar::Elements2D[ielem][ipoin] + 1 << " ";
				}
				switch (elemType)
				{
				case 3:
				{
					fileFlux << meshVar::Elements2D[ielem][0] + 1 << std::endl;
					break;
				}
				case 4:
				{
					fileFlux << meshVar::Elements2D[ielem][3] + 1 << std::endl;
					break;
				}
				default:
					break;
				}
			}
		}
		else
		{
            message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(fileName, fileLoc));
		}
	}

	namespace calcNodeValueAtBCChildFuncs
	{
		namespace patch
		{
            std::vector<double> inFlow(int element, int edgeGrp, double a, double b)
			{
				std::vector<double> U(4, 0.0), UPlus(4, 0.0), UMinus(4, 0.0);

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				//Apply weak Riemann infinite value
                UMinus[0] = (bcValues::pBCFixed[edgeGrp - 1] / (material::R*bcValues::TBCFixed[edgeGrp - 1]));
                UMinus[1] = UMinus[0] * bcValues::uBCFixed[edgeGrp - 1];
                UMinus[2] = UMinus[0] * bcValues::vBCFixed[edgeGrp - 1];
                UMinus[3] = UMinus[0] * (bcValues::TBCFixed[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBCFixed[edgeGrp - 1], 2) + pow(bcValues::vBCFixed[edgeGrp - 1], 2)));

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}

            std::vector<double> outFlow(int element, int edgeGrp, double a, double b)
			{
				std::vector<double> U(4, 0.0), UPlus(4, 0.0), UMinus(4, 0.0);

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				double uPlus(UPlus[1] / UPlus[0]), vPlus(UPlus[2] / UPlus[0]);

				//Apply PNR (2), R (1)
				int implementation(1);
				switch (implementation)
				{
				case 1: //R
				{
                    if (flowProperties::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
                        UMinus[3] = bcValues::pBCFixed[edgeGrp - 1] / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];
					}
				}
				break;
				case 2: //PNR
				{
					double TInternal(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]));
					double pInternal(0);
					pInternal = UPlus[0] * material::R*TInternal;
                    if (flowProperties::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
                        UMinus[3] = (2 * bcValues::pBCFixed[edgeGrp - 1] - pInternal) / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];
					}
				}
				break;
				default:
					break;
				}

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}
		}

		namespace wall
		{
            std::vector <double> noSlipIsoThermal(int element, int edgeGrp, double a, double b)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double> U(4, 0.0), UMinus(4, 0.0), UPlus(4, 0.0);

				//A1 approach
				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				UMinus[0] = UPlus[0];
				UMinus[1] = 0.0;
				UMinus[2] = 0.0;
                UMinus[3] = UPlus[0] * material::Cv*bcValues::TBCFixed[edgeGrp - 1];

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}

            std::vector <double> noSlipAdiabatic(int element, double a, double b)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double> U(4, 0.0), UMinus(4, 0.0), UPlus(4, 0.0);

				//A1 approach
				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				UMinus[0] = UPlus[0];
				UMinus[1] = 0.0;
				UMinus[2] = 0.0;
				UMinus[3] = UPlus[0] * material::Cv*math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]);

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}
		}

        std::vector <double> Symmetry(int element, int edge, double a, double b)
		{
			std::vector<double> U(4, 0.0), UPlus(4, 0.0), UMinus(4, 0.0);
			double nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
			for (int i = 0; i < 4; i++)
			{
				UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
			}
			UMinus[0] = UPlus[0];
			UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
			UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
			UMinus[3] = UPlus[3];
			for (int i = 0; i < 4; i++)
			{
				U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
			}
			return U;
		}
	}
}

namespace postProcessing_Surface {
    void writeVarsAtWall(std::string Loc)
    {
        int globleEdge(-1), element;
        double a, b, rhoBC;
        std::string fileName("surfaceVariables_proc");
        fileName=fileName+std::to_string(systemVar::currentProc)+".txt";
        std::ofstream FileFlux((Loc+"/"+fileName).c_str());

        if (FileFlux)
        {
            FileFlux<<message::headerFile()<<"\n"
                   <<"xCoor yCoor p T u v\n";

            for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
            {
                globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
                int bcType(auxUlti::checkBCTypeOfEdge(globleEdge));
                if (bcType==meshVar::BCTypeID::wall) //type wall
                {
                    for (int nG=0; nG<mathVar::nGauss1D+1; nG++)
                    {
                        //Tinh p tai hinh chieu vuong goc cua center xuong BCEdge
                        std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(globleEdge);
                        //Get Gauss point's coordinates
                        std::tie(a, b) = auxUlti::getGaussSurfCoor(globleEdge, element, nG);
                        rhoBC=math::pointValue(element,a,b,1,2);

                        FileFlux
                                  <<meshVar::GaussPtsOnBCEdge_x[ilocalEdge][nG]<<" "
                                  <<meshVar::GaussPtsOnBCEdge_y[ilocalEdge][nG]<<" "
                                  <<math::CalcP(SurfaceBCFields::TBc[ilocalEdge][nG],rhoBC)<<" "
                                  <<SurfaceBCFields::TBc[ilocalEdge][nG]<<" "
                                  <<SurfaceBCFields::uBc[ilocalEdge][nG]<<" "
                                  <<SurfaceBCFields::vBc[ilocalEdge][nG]<<" "
                                  <<"\n";
                        /*
                        if (auxUlti::checkTimeVaryingBCAvailable())
                            FileFlux
                                      <<meshVar::GaussPtsOnBCEdge_x[ilocalEdge][nG]<<" "
                                      <<meshVar::GaussPtsOnBCEdge_y[ilocalEdge][nG]<<" "
                                      <<math::CalcP(nonEqmSurfaceField::TBc[ilocalEdge][nG],rhoBC)<<" "
                                      <<nonEqmSurfaceField::TBc[ilocalEdge][nG]<<" "
                                      <<nonEqmSurfaceField::uBc[ilocalEdge][nG]<<" "
                                      <<nonEqmSurfaceField::vBc[ilocalEdge][nG]<<" "
                                      <<"\n";
                        else
                            FileFlux
                                      <<meshVar::GaussPtsOnBCEdge_x[ilocalEdge][nG]<<" "
                                      <<meshVar::GaussPtsOnBCEdge_y[ilocalEdge][nG]<<" "
                                      <<math::CalcP(SurfaceBCFields::TBc[ilocalEdge],rhoBC)<<" "
                                      <<SurfaceBCFields::TBc[ilocalEdge]<<" "
                                      <<SurfaceBCFields::uBc[ilocalEdge]<<" "
                                      <<SurfaceBCFields::vBc[ilocalEdge]<<" "
                                      <<"\n";
                                      */
                    }
                }
            }
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(Loc+"/"+fileName)<<" to write.\n";
        }
    }
}
