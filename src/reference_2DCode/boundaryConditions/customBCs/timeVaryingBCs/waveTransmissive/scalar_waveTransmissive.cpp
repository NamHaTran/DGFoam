#include "scalar_waveTransmissive.h"
#include "./boundaryConditions/customBCs/timeVaryingBCs/waveTransmissive/supportFuncs_waveTransmissive.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include "./boundaryConditions/bcVariables.h"
#include "./boundaryConditions/fixedValue.h"
#include <vector>
#include <iostream>

namespace waveTransmissive {
    void solveScalarEqn(int edge, int edgeGrp, int nG, std::string var)
    {
        double phi(0), Cx(0), Cy(0);
        double dPhi[2];
        int element(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(edge);
        double u=SurfaceBCFields::uBc[localEdgeId][nG], v=SurfaceBCFields::vBc[localEdgeId][nG],
                phi0(0);

        //Swith transmissive velocity between speed of sound and advective velocity
        if ((var.compare("T") == 0 && waveTransmissive::includeSoundSpeed_T[edgeGrp-1])
                || (var.compare("p") == 0 && waveTransmissive::includeSoundSpeed_p[edgeGrp-1]))
        {
            double ix(meshVar::GaussPtsOnBCEdge_unitVector_x[localEdgeId][nG]),
                    iy(meshVar::GaussPtsOnBCEdge_unitVector_y[localEdgeId][nG]);
            double a,b;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            double rhoBC(math::pointValue(element,a,b,1,2)),
                    rhouBC(math::pointValue(element,a,b,2,2)),
                    rhovBC(math::pointValue(element,a,b,3,2)),
                    rhoEBC(math::pointValue(element,a,b,4,2));

            double C(math::CalcSpeedOfSound(math::CalcTFromConsvVar(rhoBC,rhouBC,rhovBC,rhoEBC)));
            Cx=C*ix;
            Cy=C*iy;
        }

        if (var.compare("T") == 0)
            phi0=SurfaceBCFields::TBc[localEdgeId][nG];
        else if (var.compare("p") == 0)
            phi0=SurfaceBCFields::pBc[localEdgeId][nG];

        if (flowProperties::viscous)
            waveTransmissive::calcDivPhi_DG(dPhi,edge,element,nG,var);
        else
            waveTransmissive::calcDivPhi_FDM(dPhi,edge,element,nG,var);

        phi = phi0 - dt*((u+Cx)*dPhi[0]+(v+Cy)*dPhi[1]);

        if (var.compare("T") == 0)
            SurfaceBCFields::TBc[localEdgeId][nG]=phi;
        else if (var.compare("p") == 0)
            SurfaceBCFields::pBc[localEdgeId][nG]=phi;
    }

    void solveScalarEqn_implicit(int edge, int edgeGrp, int nG, std::string var)
    {
        double Cx(0), Cy(0);
        int element(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(edge);
        double u=SurfaceBCFields::uBc[localEdgeId][nG], v=SurfaceBCFields::vBc[localEdgeId][nG],
                phi0(0);

        double ix(meshVar::GaussPtsOnBCEdge_unitVector_x[localEdgeId][nG]),
                iy(meshVar::GaussPtsOnBCEdge_unitVector_y[localEdgeId][nG]),
                delta(meshVar::distanceFromGaussPtsToCentroid[localEdgeId][nG]),
                a,b;

        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
        double rhoBC(math::pointValue(element,a,b,1,2)),
                rhouBC(math::pointValue(element,a,b,2,2)),
                rhovBC(math::pointValue(element,a,b,3,2)),
                rhoEBC(math::pointValue(element,a,b,4,2));
        double rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]),
                rhoEC(rhoE[element][0]),
                phiC(0);

        //Swith transmissive velocity between speed of sound and advective velocity
        if ((var.compare("T") == 0 && waveTransmissive::includeSoundSpeed_T[edgeGrp-1])
                || (var.compare("p") == 0 && waveTransmissive::includeSoundSpeed_p[edgeGrp-1]))
        {
            double C(math::CalcSpeedOfSound(math::CalcTFromConsvVar(rhoBC,rhouBC,rhovBC,rhoEBC)));
            Cx=C*ix;
            Cy=C*iy;
        }

        if (var.compare("T") == 0)
        {
            phi0=SurfaceBCFields::TBc[localEdgeId][nG];
            phiC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);
            SurfaceBCFields::TBc[localEdgeId][nG]=(phi0+(ix*(u+Cx)+iy*(v+Cy))*dt*phiC/delta)/(1+(ix*(u+Cx)+iy*(v+Cy))*dt/delta);
        }
        else if (var.compare("p") == 0)
        {
            phi0=SurfaceBCFields::pBc[localEdgeId][nG];
            phiC=math::CalcP(math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC),rhoC);
            SurfaceBCFields::pBc[localEdgeId][nG]=(phi0+(ix*(u+Cx)+iy*(v+Cy))*dt*phiC/delta)/(1+(ix*(u+Cx)+iy*(v+Cy))*dt/delta);
        }
    }

    void calcDivPhi_DG(double *divPhi, int edge, int element, int nG, std::string var)
    {
        double a(0.0), b(0.0);
        //Tinh toa do diem Gauss dang xet
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        double rho(math::pointValue(element,a,b,1,2)),
                rhou(math::pointValue(element,a,b,2,2)),
                rhov(math::pointValue(element,a,b,3,2)),
                rhoE(math::pointValue(element,a,b,4,2)),
                drhox(math::pointAuxValue(element,a,b,1,1)),
                drhoy(math::pointAuxValue(element,a,b,1,2)),
                drhoux(math::pointAuxValue(element,a,b,2,1)),
                drhouy(math::pointAuxValue(element,a,b,2,2)),
                drhovx(math::pointAuxValue(element,a,b,3,1)),
                drhovy(math::pointAuxValue(element,a,b,3,2)),
                drhoEx(math::pointAuxValue(element,a,b,4,1)),
                drhoEy(math::pointAuxValue(element,a,b,4,2));

        double T(math::CalcTFromConsvVar(rho,rhou,rhov,rhoE));
        double mu=math::CalcVisCoef(T);

        drhox=drhox/mu;
        drhoux=drhoux/mu;
        drhovx=drhovx/mu;
        drhoEx=drhoEx/mu;
        drhox=drhox/mu;
        drhoux=drhoux/mu;
        drhovx=drhovx/mu;
        drhoEx=drhoEx/mu;

        double dux = math::calcRhouvEDeriv(drhoux, drhox, rhou, rho),
                dvx = math::calcRhouvEDeriv(drhovx, drhox, rhov, rho),
                dEx = math::calcRhouvEDeriv(drhoEx, drhox, rhoE, rho),
                dTx = math::calcTDeriv(dEx, dux, dvx, rhou/rho, rhov/rho);

        double duy = math::calcRhouvEDeriv(drhouy, drhoy, rhou, rho),
                dvy = math::calcRhouvEDeriv(drhovy, drhoy, rhov, rho),
                dEy = math::calcRhouvEDeriv(drhoEy, drhoy, rhoE, rho),
                dTy = math::calcTDeriv(dEy, duy, dvy, rhou/rho, rhov/rho);

        if (var.compare("T") == 0)
        {
            divPhi[0]=dTx;
            divPhi[1]=dTy;
        }
        else if (var.compare("p") == 0)
        {
            divPhi[0]=material::R*(dTx*rho+drhox*T);
            divPhi[1]=material::R*(dTy*rho+drhoy*T);
        }
    }

    void calcDivPhi_FDM(double *divPhi, int edge, int element, int nG, std::string var)
    {
        double a(0.0), b(0.0);
        //Tinh toa do diem Gauss dang xet
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        double TC(0),
                pC(0),
                rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]),
                rhoEC(rhoE[element][0]);

        TC = math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);
        pC = math::CalcP(TC,rhoC);

        double rhoBC(math::pointValue(element,a,b,1,2)),
                rhouBC(math::pointValue(element,a,b,2,2)),
                rhovBC(math::pointValue(element,a,b,3,2)),
                rhoEBC(math::pointValue(element,a,b,4,2));

        int localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        double delta(meshVar::distanceFromGaussPtsToCentroid[localEdgeId][nG]);
        double
                ix(meshVar::GaussPtsOnBCEdge_unitVector_x[localEdgeId][nG]),
                iy(meshVar::GaussPtsOnBCEdge_unitVector_y[localEdgeId][nG]);

        double TBC = math::CalcTFromConsvVar(rhoBC,rhouBC,rhovBC,rhoEBC),
                pBC = math::CalcP(TBC,rhoBC);

        if (var.compare("T") == 0)
        {
            divPhi[0]=(TBC-TC)*ix/delta;
            divPhi[1]=(TBC-TC)*iy/delta;
        }
        else if (var.compare("p") == 0)
        {
            divPhi[0]=(pBC-pC)*ix/delta;
            divPhi[1]=(pBC-pC)*iy/delta;
        }
    }

    /*
    void correctPhi(int edge, int edgeGrp, int nG, double &varM, double varP, std::string var)
    {
        //Apply boundary condition
        bool isStrong(BCVars::DirichletAppMethUStrong[edgeGrp -1]);
        int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        double C(0);
        if (var.compare("T") == 0)
        {
            C=SurfaceBCFields::TBc[loc][nG];
        }
        else if (var.compare("p") == 0)
        {
            C=SurfaceBCFields::pBc[loc][nG];
        }
        varM = fixedValue_scalar(varP,C,isStrong);
    }*/
}
