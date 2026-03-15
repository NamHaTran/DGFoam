#include "T_SmoluchowskyTJump.h"
#include "VarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/readBCInfor/supportReadingBCFuncs.h"
#include "./boundaryConditions/fixedValue.h"
#include <fstream>
#include <vector>
#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <math.h>
#include <iostream>

#include ".././nonEqmBCs_GenFuncs.h"

namespace SmoluchowskyTJump {
    void T_IO(int bcGrp, std::ifstream &FileFlux)
    {
        bcValues::temperatureJump=true;
        bcValues::TBcType[bcGrp - 1] = BCVars::temperatureBCId::SmoluchowskyTJump;

        std::string line, tempStr;
        std::getline(FileFlux, line);
        std::istringstream Stream1(line);
        //Read sigmaT
        Stream1>> tempStr;
        if ((tempStr.compare("sigmaT") != 0))
        {
            message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'sigmaT' in file T/group " + std::to_string(bcGrp) + ".\n");
        }
        else
        {
            Stream1>> bcValues::sigmaT;
        }

        std::getline(FileFlux, line);
        std::istringstream Stream2(line);
        //Read TWall
        Stream2 >> tempStr;
        if ((tempStr.compare("TWall") != 0))
        {
            message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'TWall' in file T/group " + std::to_string(bcGrp) + ".\n");
        }
        else
        {
            Stream2 >> bcValues::TBCFixed[bcGrp - 1];
        }

        readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethTStrong,"T");
    }

    void correctT(int edge, int edgeGrp, int nG, double &varM, double varP)
    {
        //Apply boundary condition
        bool isStrong(BCVars::DirichletAppMethTStrong[edgeGrp-1]);
        varM = fixedValue_scalar(varP,
                SurfaceBCFields::TBc[auxUlti::getAdressOfBCEdgesOnBCValsArray(edge)][nG],
                isStrong);
    }

    /*
    void correctGradT(std::vector<double> &gradM, const std::vector<double> &gradP)
    {
        gradM=gradP;
    }*/

    void calcTJump_DGTypeExplicit(int edge, int edgeGrp, int nG)
    {
        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);

        //Tinh cac gia tri can thiet-----------------------------------------------------
        //Lay cac gia tri da luu
        double aConst(0.0), dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), lambda(0.0), dEx, dEy,
                TVal(SurfaceBCFields::TBc[localEdgeId][nG]),
                muVal(0.0), nx, ny, dTx, dTy,dTn;

        double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0);
        nx = auxUlti::getNormVectorComp(element, edge, 1);
        ny = auxUlti::getNormVectorComp(element, edge, 2);

        double a(0.0), b(0.0);
        //Tinh toa do diem Gauss dang xet
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        double rhoBC(math::pointValue(element,a,b,1,2)),
                rhouBC(math::pointValue(element,a,b,2,2)),
                rhovBC(math::pointValue(element,a,b,3,2)),
                rhoEBC(math::pointValue(element,a,b,4,2)),
                dRhoBCx(math::pointAuxValue(element,a,b,1,1)),
                dRhoBCy(math::pointAuxValue(element,a,b,1,2)),
                dRhouBCx(math::pointAuxValue(element,a,b,2,1)),
                dRhouBCy(math::pointAuxValue(element,a,b,2,2)),
                dRhovBCx(math::pointAuxValue(element,a,b,3,1)),
                dRhovBCy(math::pointAuxValue(element,a,b,3,2)),
                dRhoEBCx(math::pointAuxValue(element,a,b,4,1)),
                dRhoEBCy(math::pointAuxValue(element,a,b,4,2));

        //Tinh cac thanh phan dao ham theo 2 phuong
        //dTx=math::calcDivTFromAuxVariable(dRhoEBCx,dRhoBCx,rhoEBC,rhoBC);
        //dTy=math::calcDivTFromAuxVariable(dRhoEBCy,dRhoBCy,rhoEBC,rhoBC);

        dux=math::calcRhouvEDeriv(dRhouBCx,dRhoBCx,rhouBC,rhoBC);
        duy=math::calcRhouvEDeriv(dRhouBCy,dRhoBCy,rhouBC,rhoBC);

        dvx=math::calcRhouvEDeriv(dRhovBCx,dRhoBCx,rhovBC,rhoBC);
        dvy=math::calcRhouvEDeriv(dRhovBCy,dRhoBCy,rhovBC,rhoBC);

        dEx=math::calcRhouvEDeriv(dRhoEBCx,dRhoBCx,rhoEBC,rhoBC);
        dEy=math::calcRhouvEDeriv(dRhoEBCy,dRhoBCy,rhoEBC,rhoBC);

        dTx=math::calcTDeriv(dEx,dux,dvx,rhouBC/rhoBC,rhovBC/rhoBC);
        dTy=math::calcTDeriv(dEy,duy,dvy,rhouBC/rhoBC,rhovBC/rhoBC);

        //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
        dTn=dTx*nx+dTy*ny;
        muVal=math::CalcVisCoef(TVal);

        lambda=math::calcMeanFreePath(muVal,rhoBC,TVal)/muVal;
        aConst=(2-bcValues::sigmaT)*2*material::gamma*lambda/(bcValues::sigmaT*material::Pr*(material::gamma+1));

        TJump=TWall-aConst*dTn;
        if (TJump<0)
        {
            TJump=TVal;
        }
        //Update to surfaceBCfields
        SurfaceBCFields::TBc[localEdgeId][nG] = TJump;
    }

    void calcTJump_FDMTypeImplicit(int edge, int edgeGrp, int nG)
    {
        //Dieu kien bien Maxwell-Smoluchowski

        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);
        //Tinh cac gia tri can thiet
        double TC, rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]),
                rhoEC(rhoE[element][0]),
                TVal(SurfaceBCFields::TBc[localEdgeId][nG]),
                delta(meshVar::distanceFromGaussPtsToCentroid[localEdgeId][nG]);

        double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0), a, b;

        double
                nx = auxUlti::getNormVectorComp(element, edge, 1),
                ny = auxUlti::getNormVectorComp(element, edge, 2);

        //Lay vector don vi (C)->(Gauss)
        double
                ix(meshVar::GaussPtsOnBCEdge_unitVector_x[localEdgeId][nG]),
                iy(meshVar::GaussPtsOnBCEdge_unitVector_y[localEdgeId][nG]);

        //Tinh dot product i.n
        double iDotn(nx*ix+ny*iy);

        //Tinh toa do diem Gauss dang xet
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        double rhoBC(math::pointValue(element,a,b,1,2));

        //Gia tri trung binh (tai centroid)
        TC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);

        //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
        /* Note: term dT/dn, normal gradient component of T, is approximated by backward finite difference method
         *                                  dT/dn = ((TJump - TC)/delta)*(i dot n)
         * with
         *      TJump: temperature at wall (face of cell)
         *      TC: temperature at cell centroid
         *      delta: distance from centroid to surface Gauss point
         *      i dot n: dot product of i and n
         */
        double AT(material::viscosityCoeff::Sutherland::As*pow(3.1416/(2*material::R),0.5)/rhoBC),
                BT(2*material::gamma*(2-bcValues::sigmaT)/(bcValues::sigmaT*material::Pr*(material::gamma+1)));
        double ATT(delta+AT*BT*iDotn),
                BTT(material::viscosityCoeff::Sutherland::Ts*delta-TWall*delta-AT*BT*TC*iDotn),
                CTT(-TWall*material::viscosityCoeff::Sutherland::Ts*delta), T1, T2;

        bool realSolution(false);

        std::tie(realSolution, T1, T2) = math::solvePolynomialsEq::polynominal2d(ATT,BTT,CTT);
        if (realSolution)
        {
            if (T1>0)
            {
                TJump=T1;
            }
            else
            {
                TJump=T2;
            }
        }
        else
        {
            TJump=TVal;
        }
        //Update to surfaceBCfields
        SurfaceBCFields::TBc[localEdgeId][nG] = TJump;
    }


    void calcTJump_FDMTypeSemiImplicit(int edge, int edgeGrp, int nG)
    {
        //Dieu kien bien Maxwell-Smoluchowski

        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);
        //Tinh cac gia tri can thiet
        double TC, rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]),
                rhoEC(rhoE[element][0]),
                TVal(SurfaceBCFields::TBc[localEdgeId][nG]),
                delta(meshVar::distanceFromGaussPtsToCentroid[localEdgeId][nG]);

        double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(SurfaceBCFields::TBc[localEdgeId][nG]), a, b;

        double
                nx = auxUlti::getNormVectorComp(element, edge, 1),
                ny = auxUlti::getNormVectorComp(element, edge, 2);

        //Lay vector don vi (C)->(Gauss)
        double
                ix(meshVar::GaussPtsOnBCEdge_unitVector_x[localEdgeId][nG]),
                iy(meshVar::GaussPtsOnBCEdge_unitVector_y[localEdgeId][nG]);

        //Tinh dot product i.n
        double iDotn(nx*ix+ny*iy);

        //Tinh toa do diem Gauss dang xet
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        double rhoBC(math::pointValue(element,a,b,1,2));

        //Gia tri trung binh (tai centroid)
        TC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);

        //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
        /* Note: term dT/dn, normal gradient component of T, is approximated by backward finite difference method
         *                                  dT/dn = ((TJump - TC)/delta)*(i dot n)
         * with
         *      TJump: temperature at wall (face of cell)
         *      TC: temperature at cell centroid
         *      delta: distance from centroid to surface Gauss point
         *      i dot n: dot product of i and n
         */
        double lambda(0.0), A(0.0);

        lambda=math::calcMeanFreePath(
                    math::CalcVisCoef(TJump),
                    rhoBC,
                    TJump
                    );
        A=(2-bcValues::sigmaT)*2*material::gamma*lambda/(bcValues::sigmaT*(material::gamma+1)*material::Pr)*iDotn/delta;

        TJump=(TWall+A*TC)/(1+A);

        if (TJump<0)
        {
            TJump=TVal;
        }
        //Update to surfaceBCfields
        SurfaceBCFields::TBc[localEdgeId][nG] = TJump;
    }
}
