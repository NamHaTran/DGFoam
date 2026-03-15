#include "u_MaxwellSlip.h"
#include "VarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/readBCInfor/supportReadingBCFuncs.h"
#include <fstream>
#include "./boundaryConditions/fixedValue.h"
#include <vector>
#include "DGMath.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <math.h>
#include <iostream>

#include ".././nonEqmBCs_GenFuncs.h"

//Extended Navier-Stokes-Fourier model
#include "./extNSFEqns/FranzDurst/DurstModel.h"

namespace MaxwellSlip {
    void u_IO(int bcGrp, std::ifstream &FileFlux)
    {
        bcValues::slipBCFlag=true;
        bcValues::UBcType[bcGrp - 1] = BCVars::velocityBCId::MaxwellSlip;

        std::string line, tempStr;
        std::getline(FileFlux, line);
        std::istringstream fixedUStream1(line);
        //Read sigmaU
        fixedUStream1>>tempStr;
        if ((tempStr.compare("sigmaU") != 0))
        {
            message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'sigmaU' in file U/group " + std::to_string(bcGrp) + ".\n");
        }
        else
        {
            fixedUStream1>>bcValues::sigmaU;
        }

        std::getline(FileFlux, line);
        std::istringstream fixedUStream2(line);
        //Read UWall
        fixedUStream2 >> tempStr;
        if ((tempStr.compare("uWall") != 0))
        {
            message::writeLog(systemVar::pwd, systemVar::caseName, "Cannot find key word 'uWall' in file U/group " + std::to_string(bcGrp) + ".\n");
        }
        else
        {
            fixedUStream2 >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
        }

        //read U/gradU application method
        readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong,"U");
    }

    void correctU(int edge, int edgeGrp, int nG, std::vector<double> &varM, std::vector<double> varP)
    {
        //Apply boundary condition
        bool isStrong(BCVars::DirichletAppMethUStrong[edgeGrp -1]);
        std::vector<double> UBC{
                    SurfaceBCFields::uBc[auxUlti::getAdressOfBCEdgesOnBCValsArray(edge)][nG],    //uSlip
                    SurfaceBCFields::vBc[auxUlti::getAdressOfBCEdgesOnBCValsArray(edge)][nG]};   //vSlip
        fixedValue_vector(varM,varP,UBC,isStrong);
    }

    /*
    void correctGradU(std::vector<double> &gradM, const std::vector<double> &gradP)
    {
        gradM=gradP;
    }*/

    void calcUSlip_DGTypeExplicit(int edge, int edgeGrp, int nG)
    {
        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);

        //Tinh cac gia tri can thiet-----------------------------------------------------
        //Lay cac gia tri da luu
        double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), lambda(0.0), dEx, dEy,
                muVal(0.0), nx, ny, dTx, dTy;

        double TJump(SurfaceBCFields::TBc[localEdgeId][nG]),
        uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]);
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

        //Maxwell slip velocity--------------------------------------------------------
        //Tinh lai gia tri muVal va lambda tu TJump
        muVal=math::CalcVisCoef(TJump);
        lambda=math::calcMeanFreePath(muVal,rhoBC,TJump);

        //Trasnformation tensor components
        double S11(1-nx*nx),
                S12(-nx*ny),
                S21(-nx*ny),
                S22(1-ny*ny);

        //- Term div_n(S.u)
        double div_n_u1(dux*nx+dvx*ny), div_n_u2(duy*nx+dvy*ny);
        std::vector<double> div_n_Sdotu={
            S11*div_n_u1 + S12*div_n_u2,
            S21*div_n_u1 + S22*div_n_u2
        };

        //- Term S.(n.PI_mc)
        double PImc_n1(muVal*(nx*(dux-(2/3)*(dux+dvy))+ny*duy)), PImc_n2(muVal*(nx*dvx+ny*(dvy-(2/3)*(dux+dvy))));
        std::vector<double> SdotndotPImc={
            S11*PImc_n1 + S12*PImc_n2,
            S21*PImc_n1 + S22*PImc_n2
        };

        //- Term S.div(T)
        std::vector<double> SdotDivT={
            S11*dTx + S12*dTy,
            S21*dTx + S22*dTy
        };

        //Giai slip velocity
        double c=(2-bcValues::sigmaU)/bcValues::sigmaU;
        double
                uSlip=uWall-c*lambda*SdotndotPImc[0]/muVal-(3/4)*SdotDivT[0]*muVal/(rhoBC*TJump)-c*lambda*div_n_Sdotu[0],
                vSlip=vWall-c*lambda*SdotndotPImc[1]/muVal-(3/4)*SdotDivT[1]*muVal/(rhoBC*TJump)-c*lambda*div_n_Sdotu[1];

        //Calculate diffusive velocity
        double uD(0.0), vD(0.0);
        if (flowProperties::massDiffusion)
        {
            std::vector<double> U {
                math::pointValue(element,a,b,1,2),
                math::pointValue(element,a,b,2,2),
                math::pointValue(element,a,b,3,2),
                math::pointValue(element,a,b,4,2)},
            dUx {
                math::pointAuxValue(element,a,b,1,1),
                math::pointAuxValue(element,a,b,2,1),
                math::pointAuxValue(element,a,b,3,1),
                math::pointAuxValue(element,a,b,4,1)},
            dUy {
                math::pointAuxValue(element,a,b,1,2),
                math::pointAuxValue(element,a,b,2,2),
                math::pointAuxValue(element,a,b,3,2),
                math::pointAuxValue(element,a,b,4,2)};

            uD=extNSF_Durst::calcDiffVelocity(U,dUx);
            vD=extNSF_Durst::calcDiffVelocity(U,dUy);
        }

        //Update arrays
        SurfaceBCFields::uBc[localEdgeId][nG] = uSlip + uD;
        SurfaceBCFields::vBc[localEdgeId][nG] = vSlip + vD;
    }

    void calcUSlip_FDMTypeImplicit(int edge, int edgeGrp, int nG)
    {
        //Dieu kien bien Maxwell-Smoluchowski

        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);
        //Tinh cac gia tri can thiet
        double lambda(0.0),
                TC, uC, vC,
                rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]),
                rhoEC(rhoE[element][0]),
                delta(meshVar::distanceFromGaussPtsToCentroid[localEdgeId][nG]),
                muVal(0.0), nx, ny, dTx, dTy,dTn;

        double TJump(0.0), uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]), a, b;

        //Tinh toa do diem Gauss dang xet
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        double rhoBC(math::pointValue(element,a,b,1,2));
        nx = auxUlti::getNormVectorComp(element, edge, 1);
        ny = auxUlti::getNormVectorComp(element, edge, 2);

        //Lay vector don vi (C)->(Gauss)
        double
                ix(meshVar::GaussPtsOnBCEdge_unitVector_x[localEdgeId][nG]),
                iy(meshVar::GaussPtsOnBCEdge_unitVector_y[localEdgeId][nG]);

        //Tinh dot product i.n
        double iDotn(nx*ix+ny*iy);

        //Gia tri trung binh (tai centroid)
        TC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);
        uC=rhouC/rhoC;
        vC=rhovC/rhoC;

        //Lay gia tri TJump trong array
        TJump=SurfaceBCFields::TBc[localEdgeId][nG];

        //Gia tri normal gradient (by following Finite Difference Method)
        dTn=(TJump-TC)*iDotn/delta;

        //Tinh cac thanh phan dao ham theo 2 phuong
        dTx=dTn*nx;
        dTy=dTn*ny;

        muVal=math::CalcVisCoef(TJump);
        lambda=math::calcMeanFreePath(muVal,rhoBC,TJump);

        //2. Maxwell slip velocity--------------------------------------------------------
        /* Note: term dU/dn, normal gradient component of U(u, v), is approximated by backward finite difference method
         *                                  dU/dn = ((USlip - UC)/delta)*(i dot n)
         * with
         *      USlip: velocity at wall (face of cell)
         *      UC: velocity at cell centroid
         *      delta: distance from centroid to surface Gauss point
         *      i dot n: dot product of i and n
         */

        uC=uC/delta;
        vC=vC/delta;

        //Transformation tensor components
        double S11(1-nx*nx),
                S12(-nx*ny),
                S21(-nx*ny),
                S22(1-ny*ny);

        double A1=(2-bcValues::sigmaU)*lambda/bcValues::sigmaU;
        //- Term (-3/4)*mu*S.div(T)/(rho*T)
        std::vector<double> A2={
            -3*muVal*(S11*dTx + S12*dTy)/(4*rhoBC*TJump),
            -3*muVal*(S21*dTx + S22*dTy)/(4*rhoBC*TJump)
        };

        //He so B
        double B1(-(uC*nx*nx+vC*nx*ny)*iDotn),
                B2(-(uC*nx*ny+vC*ny*ny)*iDotn),
                B3(nx*nx*iDotn/delta),
                B4(nx*ny*iDotn/delta),
                B5(ny*ny*iDotn/delta);

        //He so C
        double C1(S11*B3+S12*B4),
                C2(S11*B4+S12*B5),
                C3(S11*B1+S12*B2),
                C4(S21*B3+S22*B4),
                C5(S21*B4+S22*B5),
                C6(S21*B1+S22*B2);

        //He so D
        double D1(iDotn*(muVal*(nx*nx/(3*delta)+ny*ny/delta))),
                D2(iDotn*(muVal*(-2*nx*ny/(3*delta)))),
                D3(iDotn*(muVal*(-uC*ny*ny-uC*nx*nx/3+2*nx*ny*vC/3))),
                D4(iDotn*(muVal*(-2*nx*ny/(3*delta)))),
                D5(iDotn*(muVal*(ny*ny/(3*delta)+nx*nx/delta))),
                D6(iDotn*(muVal*(-vC*nx*nx-vC*ny*ny/3+2*nx*ny*uC/3)));

        //He so E
        double E1(S11*D1+S12*D4),
                E2(S11*D2+S12*D5),
                E3(S11*D3+S12*D6),
                E4(S21*D1+S22*D4),
                E5(S21*D2+S22*D5),
                E6(S21*D3+S22*D6);

        //He so F
        double F1(1+A1*C1+A1*E1/muVal),
                F2(A1*C2+A1*E2/muVal),
                F3(uWall-(A1/muVal)*E3-A1*C3+A2[0]),
                F4(A1*C4+A1*E4/muVal),
                F5(1+A1*C5+A1*E5/muVal),
                F6(vWall-(A1/muVal)*E6-A1*C6+A2[1]);

        //Giai uSlip, vSlip
        double uSlip, vSlip;
        vSlip=(F6-F4*F3/F1)/((F5-F4*F2/F1));
        uSlip=(F3-F2*vSlip)/(F1);

        //Calculate diffusive velocity
        double uD(0.0), vD(0.0);
        if (flowProperties::massDiffusion)
        {
            //uD=extNSF_Durst::uMassWall[localEdgeId];
            //vD=extNSF_Durst::vMassWall[localEdgeId];
        }

        //Update
        /* Day 21/04/2022:
         * Khong tinh them thanh phan van toc diffusive vao trong van toc truot vi:
         * - Bien U trong code la U_convective, nen U_slip, se duoc dung de lam dk bien, cung phai la
         * U_convective.
        */
        SurfaceBCFields::uBc[localEdgeId][nG] = uSlip + uD;
        SurfaceBCFields::vBc[localEdgeId][nG] = vSlip + vD;
    }

    double calcGradientUsingFDM(double phiC, double phiBC, double delta, double iComponent)
    {
        return (((phiBC - phiC)/delta)*iComponent);
    }
}
