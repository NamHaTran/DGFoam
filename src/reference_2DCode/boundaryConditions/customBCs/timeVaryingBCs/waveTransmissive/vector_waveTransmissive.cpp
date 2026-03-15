#include "vector_waveTransmissive.h"
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
    void solveVectorEqn(int edge, int edgeGrp, int nG)
    {
        double dU[4];
        double Cx(0), Cy(0);
        int element(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(edge);

        //Swith transmissive velocity between speed of sound and advective velocity
        if (waveTransmissive::includeSoundSpeed_u[edgeGrp-1])
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

        if (flowProperties::viscous)
            waveTransmissive::calcDivU_DG(dU,edge,element,nG);
        else
            waveTransmissive::calcDivU_FDM(dU,edge,element,nG);

        double ux(dU[0]), uy(dU[1]), vx(dU[2]), vy(dU[3]);
        double U0[2] = {SurfaceBCFields::uBc[localEdgeId][nG]-dt*Cx*ux-dt*Cy*uy, SurfaceBCFields::vBc[localEdgeId][nG]-dt*Cx*vx-dt*Cy*vy};

        double uBC(0), vBC(0);
        dU[0] = dt*ux + 1;
        dU[1] = dt*uy;
        dU[2] = dt*vx;
        dU[3] = dt*vy + 1;

        std::tie(uBC,vBC) = math::solveSys2Eqs(dU,U0);

        SurfaceBCFields::uBc[localEdgeId][nG] = uBC;
        SurfaceBCFields::vBc[localEdgeId][nG] = vBC;
    }

    void calcDivU_DG(double *divU, int edge, int element, int nG)
    {
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
                dRhovBCy(math::pointAuxValue(element,a,b,3,2));

        double muBC=math::CalcVisCoef(math::CalcTFromConsvVar(rhoBC,rhouBC,rhovBC,rhoEBC));

        divU[0]=math::calcRhouvEDeriv(dRhouBCx,dRhoBCx,rhouBC,rhoBC)/muBC; //dux
        divU[1]=math::calcRhouvEDeriv(dRhouBCy,dRhoBCy,rhouBC,rhoBC)/muBC; //duy

        divU[2]=math::calcRhouvEDeriv(dRhovBCx,dRhoBCx,rhovBC,rhoBC)/muBC; //dvx
        divU[3]=math::calcRhouvEDeriv(dRhovBCy,dRhoBCy,rhovBC,rhoBC)/muBC; //dvy
    }

    void calcDivU_FDM(double *divU, int edge, int element, int nG)
    {
        double a(0.0), b(0.0);
        //Tinh toa do diem Gauss dang xet
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        double uC, vC,
                rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]);

        uC=rhouC/rhoC;
        vC=rhovC/rhoC;

        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
        double rhoBC(math::pointValue(element,a,b,1,2)),
                rhouBC(math::pointValue(element,a,b,2,2)),
                rhovBC(math::pointValue(element,a,b,3,2));

        int localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        double delta(meshVar::distanceFromGaussPtsToCentroid[localEdgeId][nG]);
        double
                ix(meshVar::GaussPtsOnBCEdge_unitVector_x[localEdgeId][nG]),
                iy(meshVar::GaussPtsOnBCEdge_unitVector_y[localEdgeId][nG]);

        double uBC=rhouBC/rhoBC,vBC=rhovBC/rhoBC;

        divU[0]=(uBC-uC)*ix/delta; //dux
        divU[1]=(uBC-uC)*iy/delta; //duy

        divU[2]=(vBC-vC)*ix/delta; //dvx
        divU[3]=(vBC-vC)*iy/delta; //dvy
    }

    /*
    void correctU(int edge, int edgeGrp, int nG, std::vector<double> &varM, std::vector<double> varP)
    {
        //Apply boundary condition
        bool isStrong(BCVars::DirichletAppMethUStrong[edgeGrp -1]);
        std::vector<double> UBC{
                    SurfaceBCFields::uBc[auxUlti::getAdressOfBCEdgesOnBCValsArray(edge)][nG],
                    SurfaceBCFields::vBc[auxUlti::getAdressOfBCEdgesOnBCValsArray(edge)][nG]};
        fixedValue_vector(varM,varP,UBC,isStrong);
    }*/
}
