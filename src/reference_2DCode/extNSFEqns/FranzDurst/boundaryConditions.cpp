#include "./extNSFEqns/FranzDurst/boundaryConditions.h"
#include <tuple>
#include <vector>
#include "DGAuxUltilitiesLib.h"
#include "dynamicVarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "DGMath.h"
#include "./extNSFEqns/FranzDurst/DurstModel.h"

namespace bcForExtNSF_Durst {
double tempdTx(0.0), tempdTy(0.0);

    void checkConditionToAddDiffTerms(int edge)
    {
        int edgeGrp(auxUlti::getGrpOfEdge(edge));
        int UType(bcValues::UBcType[edgeGrp - 1]);

        //Neu khong cho self diffusion tai wall
        if (!extNSF_Durst::diffusionAtWall)
        {
            //Neu edge dang set la wall
            if (UType==BCVars::velocityBCId::noSlip
                    || UType==BCVars::velocityBCId::movingWall
                    || UType==BCVars::velocityBCId::slipWall)
                //Thi remove diff terms
                extNSF_Durst::needToRemoveDiffTerm=true;
            else
                extNSF_Durst::needToRemoveDiffTerm=false;
        }
    }

    void checkConditionToDropNormSelfDiffTerm(int edge)
    {
        //Neu cho self diffusion tai wall
        if (extNSF_Durst::diffusionAtWall)
        {
            //Neu edge dang set la wall thi remove diff terms
            if (auxUlti::checkBCTypeOfEdge(edge) == meshVar::BCTypeID::wall)
                extNSF_Durst::dropNormSelfDiffTerm=true;
            else
                extNSF_Durst::dropNormSelfDiffTerm=false;
        }
    }

    void resetNeedToRemoveDiffTermFlag()
    {
        //Mac du khong can phai tao function nhung van tao de tranh viec quen reset
        extNSF_Durst::needToRemoveDiffTerm=false;
    }

    void resetNeedToDropNormSelfDiffTermFlag()
    {
        //Mac du khong can phai tao function nhung van tao de tranh viec quen reset
        extNSF_Durst::dropNormSelfDiffTerm=false;
    }

    void correctViscousTerms(std::vector<std::vector<double>> &diffTerms, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy, std::vector<double> &n)
    {
        std::vector<std::vector<double>> selfDiffTensor(2, std::vector<double>(4, 0.0));
        selfDiffTensor=bcForExtNSF_Durst::calcSelfDiffusionTensor(U, dUx, dUy, n);

        std::tie(diffTerms[0][0], diffTerms[0][1], diffTerms[0][2], diffTerms[0][3]) = extNSF_Durst::calcSelfDiffusionTerms(selfDiffTensor,U[1]/U[0], U[2]/U[0], 1);
        std::tie(diffTerms[1][0], diffTerms[1][1], diffTerms[1][2], diffTerms[1][3]) = extNSF_Durst::calcSelfDiffusionTerms(selfDiffTensor,U[1]/U[0], U[2]/U[0], 2);
    }

    std::vector<std::vector<double>> calcSelfDiffusionTensor(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy, std::vector<double> &n)
    {
        /*Output matrix has form:
        [mDx    tauXx		tauXy		Qx]
        [mDy    tauYx		tauYy		Qy]
        */
        std::vector<std::vector<double>> OutputMatrix(2, std::vector<double>(4, 0.0));

        double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), dEx(0.0), dEy(0.0),
                rhouVal(0.0), rhovVal, rhoVal(0.0), rhoEVal(0.0), TVal;
        double drhox(0.0), drhoy(0.0),
            drhoux(0.0), drhouy(0.0),
            drhovx(0.0), drhovy(0.0),
            drhoEx(0.0), drhoEy(0.0),
            dTx(0.0), dTy(0.0);
        double mDx(0.0), mDy(0.0);

        double uVal(0.0), vVal(0.0), eVal(0.0);
        int index(0);

        rhoVal = U[0];
        rhouVal = U[1];
        rhovVal = U[2];
        rhoEVal = U[3];

        TVal = math::CalcTFromConsvVar(rhoVal,rhouVal,rhovVal,rhoEVal);
        eVal = material::Cv*TVal;

        uVal = rhouVal / rhoVal;
        vVal = rhovVal / rhoVal;

        drhox = dUx[0];
        drhoy = dUy[0];

        drhoux = dUx[1];
        drhouy = dUy[1];

        drhovx = dUx[2];
        drhovy = dUy[2];

        drhoEx = dUx[3];
        drhoEy = dUy[3];

        dux = math::calcRhouvEDeriv(drhoux, drhox, rhouVal, rhoVal);
        duy = math::calcRhouvEDeriv(drhouy, drhoy, rhouVal, rhoVal);

        dvx = math::calcRhouvEDeriv(drhovx, drhox, rhovVal, rhoVal);
        dvy = math::calcRhouvEDeriv(drhovy, drhoy, rhovVal, rhoVal);

        dEx = math::calcRhouvEDeriv(drhoEx, drhox, rhoEVal, rhoVal);
        dEy = math::calcRhouvEDeriv(drhoEy, drhoy, rhoEVal, rhoVal);

        dTx = math::calcTDeriv(dEx, dux, dvx, uVal, vVal);
        dTy = math::calcTDeriv(dEy, duy, dvy, uVal, vVal);

        //Remove Sorret Term at Wall
        if (extNSF_Durst::dropNormSelfDiffTerm)
        {
            dTx=0;
            dTy=0;
        }

        /*selfDiffusionTensor:
        [mDx    tauXx		tauXy		Qx]
        [mDy    tauYx		tauYy		Qy]

        Details:
        mDx             (selfDiffusionTensor[0][0])
        mDy             (selfDiffusionTensor[1][0])
        tauXy           (selfDiffusionTensor[0][2])
        tauXx           (selfDiffusionTensor[0][1])
        tauYy           (selfDiffusionTensor[1][2])
        Qx              (selfDiffusionTensor[0][3])
        Qy              (selfDiffusionTensor[1][3])
        */
        //Luu y cac dao ham deu da nhan mu
        //mDx = -extNSF_Durst::Dm*(drhox/rhoVal + dTx/(TVal));
        //mDy = -extNSF_Durst::Dm*(drhoy/rhoVal + dTy/(TVal));
        mDx = extNSF_Durst::calcSelfDiffFlux(rhoVal,TVal,drhox,dTx);
        mDy = extNSF_Durst::calcSelfDiffFlux(rhoVal,TVal,drhoy,dTy);

        if (extNSF_Durst::dropNormSelfDiffTerm)
        {
            bcForExtNSF_Durst::removeNormTermOfVector(mDx,mDy,n,mDx,mDy,false);
        }

        OutputMatrix[0][0]=mDx;
        OutputMatrix[1][0]=mDy;

        /*calculate stresses*/
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                index = 10 * (i + 1) + (j + 1);
                if (index == 11)  //x normal stress (tau_xx)
                {
                    OutputMatrix[i][j+1] = extNSF_Durst::calcDiffusionStressComponent(index, mDx*uVal, mDy*vVal);
                }
                else if (index == 22)  //y normal stress (tau_yy)
                {
                    OutputMatrix[i][j+1] = extNSF_Durst::calcDiffusionStressComponent(index, mDy*vVal, mDx*uVal);
                }
                else  //shear stress (tau_xy)
                {
                    OutputMatrix[i][j+1] = extNSF_Durst::calcDiffusionStressComponent(index, mDx*vVal, mDy*uVal);
                }
            }
        }

        /*Sel-diffusion molecular transport of heat*/
        OutputMatrix[0][3] = mDx*eVal;
        OutputMatrix[1][3] = mDy*eVal;

        return OutputMatrix;
    }

    void removeNormTermOfVector(double &dx, double &dy, const std::vector<double> &n, const double vectorDirX, const double vectorDirY, bool removeOnlyWhenVectorDirIsSameDirectionWithn)
    {
        double tempX, tempY, nx(n[0]), ny(n[1]);

        if (!removeOnlyWhenVectorDirIsSameDirectionWithn || (removeOnlyWhenVectorDirIsSameDirectionWithn && (vectorDirX*nx+vectorDirY*ny)>0))
        {
            tempX = dx - 2*(dx * nx + dy * ny)*nx;
            tempY = dy - 2*(dx * nx + dy * ny)*ny;

            dx=tempX;
            dy=tempY;
        }
    }
}
