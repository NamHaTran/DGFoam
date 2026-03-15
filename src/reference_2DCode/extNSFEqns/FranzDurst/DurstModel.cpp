#include "./extNSFEqns/FranzDurst/DurstModel.h"
#include <math.h>
#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <tuple>


/* GHI CHU:
 * Cac ham trong source code chinh da modified de include model Durst la:
 * process::NSFEq::calcLocalViscousFlux;                    ------> tinh local viscous term (dung cho surface)
 * process::NSFEq::calcVolumeIntegralTerms;                 ------> them term tinh Tich phan cua Kinetic Energy Flux do self-diffusion
 * NSFEqBCsForNormalBoundary trong ham NSFEqBCsImplement    ------> setup flag de xac dinh co turn of self-diffusion tai wall khong
 * process::NSFEq::calcSurfaceFlux                          ------> setup flag de xac dinh co turn of self-diffusion tai wall khong
 * process::NSFEq::calcVolumeIntegralTerms                  ------> setup flag de xac dinh co turn of self-diffusion tai wall khong
 * message::showCaseInformations                            ------> Show Durst mode settings
 * IO::loadSettingFiles::massDiffSettings                   ------> Load settings
 * symmetryBC_NSFEqs                                        ------> setup flag de xac dinh BC dang tinh la symmetry. Flag nay de turn on ham correct flux tai symmetry
 *
 *
 * Keyword trong code: DURST MODEL
*/
namespace extNSF_Durst {
    /*Variables*/
    bool enable(false),
        diffusionAtWall(false),
        needToRemoveDiffTerm(false),

        massDiffModel_ChapmanEnskog(false),
        massDiffModel_constant(false),

        dropNormSelfDiffTerm(false),

        isSymmetry(false);

    double Dm(1.32),
        blending(0.5),
        realDm(1.32);

    double *uMassWall = new double[1];
    double *vMassWall = new double[1];

    void resizeArrays()
    {
        extNSF_Durst::uMassWall = new double [meshVar::numBCEdges];
        auxUlti::initialize1DArray(extNSF_Durst::uMassWall, meshVar::numBCEdges, 0.0);

        extNSF_Durst::vMassWall = new double [meshVar::numBCEdges];
        auxUlti::initialize1DArray(extNSF_Durst::vMassWall, meshVar::numBCEdges, 0.0);
    }

    void applyBlendingFactorToDm()
    {
        extNSF_Durst::realDm=(extNSF_Durst::Dm - 1.0)*extNSF_Durst::blending + 1.0;
    }

    void correctViscousTerms(std::vector<std::vector<double>> &diffTerms, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
    {
        std::vector<std::vector<double>> selfDiffTensor(2, std::vector<double>(4, 0.0));
        selfDiffTensor=extNSF_Durst::calcSelfDiffusionTensor(U, dUx, dUy);

        std::tie(diffTerms[0][0], diffTerms[0][1], diffTerms[0][2], diffTerms[0][3]) = extNSF_Durst::calcSelfDiffusionTerms(selfDiffTensor,U[1]/U[0], U[2]/U[0], 1);
        std::tie(diffTerms[1][0], diffTerms[1][1], diffTerms[1][2], diffTerms[1][3]) = extNSF_Durst::calcSelfDiffusionTerms(selfDiffTensor,U[1]/U[0], U[2]/U[0], 2);
    }

    std::tuple<double, double, double, double> calcSelfDiffusionTerms(std::vector< std::vector<double> > &selfDiffusionTensor, double uVal, double vVal, int dir)
    {
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

        double mDx(selfDiffusionTensor[0][0]),
                mDy(selfDiffusionTensor[1][0]),
                tauXy(selfDiffusionTensor[0][2]),
                tauXx(selfDiffusionTensor[0][1]),
                tauYy(selfDiffusionTensor[1][2]),
                Qx(selfDiffusionTensor[0][3]),
                Qy(selfDiffusionTensor[1][3]);
        double viscTerm1(0.0), viscTerm2(0.0), viscTerm3(0.0), viscTerm4(0.0);

        if (dir==1)
        {
            /*1. Ox direction*/
            viscTerm1 = mDx;
            viscTerm2 = tauXx;
            viscTerm3 = tauXy;
            viscTerm4 = tauXx * uVal + tauXy * vVal + Qx - mDx*0.5*(uVal*uVal + vVal*vVal);
        }
        else if (dir==2)
        {
            /*2. Oy direction*/
            viscTerm1 = mDy;
            viscTerm2 = tauXy;
            viscTerm3 = tauYy;
            viscTerm4 = tauXy * uVal + tauYy * vVal + Qy - mDy*0.5*(uVal*uVal + vVal*vVal);
        }
        return std::make_tuple(viscTerm1, viscTerm2, viscTerm3, viscTerm4);
    }

    std::vector<std::vector<double>> calcSelfDiffusionTensor(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
    {
        /*Output matrix has form:
        [mDx    tauXx		tauXy		Qx]
        [mDy    tauYx		tauYy		Qy]
        */
        std::vector<std::vector<double>> OutputMatrix(2, std::vector<double>(4, 0.0));

        double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), dEx(0.0), dEy(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0), rhoEVal(0.0), TVal;
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

    double calcDiffusionStressComponent(int index, double fstTerm, double sndTerm)
    {
        /*Formula of stress
        tau_xx = (4/3)*mDx*ux - (2/3)mDy*uy
        tau_yy = (4/3)*mDy*uy - (2/3)mDx*ux
        tau_xy = mDx*uy + mDy*ux

        input variables:
        index			fstTerm         sndTerm
        11 (tau_xx)		mDx*ux			mDy*Uy
        22 (tau_yy)		mDy*Uy			mDx*ux
        12 (tau_xy)		mDx*uy			mDy*ux

        Dau - trong cac term deu da co trong mDx mDy
        */
        double tau(0.0);
        if (index == 11 || index == 22)
        {
            tau = (4.0 / 3.0)*fstTerm - (2.0 / 3.0)*sndTerm;
        }
        else if (index == 12 || index == 21)
        {
            tau = fstTerm + sndTerm;
        }
        return tau;
    }

    void correctEnergyEqnVolIntTerm(int element, std::vector<double> &EnergyEqnVolIntTerm)
    {
        double a(0.0), b(0.0);
        double Int(0.0), EcOrder(0.0), dBx(0.0), dBy(0.0),
                rhoVal, rhouVal, rhovVal, rhoEVal, uVal, vVal, TVal,
                drhox, drhoy, drhoux, drhouy, drhovx, drhovy, drhoEx, drhoEy,
                dux, duy, dvx, dvy, dEx, dEy, dTx, dTy;
        int elemType(auxUlti::checkType(element));

        std::vector<double> dUx(4, 0.0),
                        dUy(4, 0.0);

        std::vector<std::vector<double>>
                mDx(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
                mDy(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));

        //Calculate Self-Diffusion Flux at all Gauss Points
        for (int na = 0; na <= mathVar::nGauss2D; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
            {
                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                std::tie(a, b) = auxUlti::getGaussCoor(na, nb);

                rhoVal=volumeFields::rhoVolGauss[element][nanb];
                rhouVal=volumeFields::rhouVolGauss[element][nanb];
                rhovVal=volumeFields::rhovVolGauss[element][nanb];
                rhoEVal=volumeFields::rhoEVolGauss[element][nanb];

                uVal=rhouVal/rhoVal;
                vVal=rhovVal/rhoVal;
                TVal=math::CalcTFromConsvVar(rhoVal,rhouVal,rhovVal,rhoEVal);

                //Tinh dao ham
                dUx=math::pointSVars(0,element,a,b,1,1);
                dUy=math::pointSVars(0,element,a,b,2,1);

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

                //Term mD
                //mDx[na][nb] = -extNSF_Durst::Dm*(drhox/rhoVal + dTx/(TVal));
                //mDy[na][nb] = -extNSF_Durst::Dm*(drhoy/rhoVal + dTy/(TVal));
                mDx[na][nb] = extNSF_Durst::calcSelfDiffFlux(rhoVal,TVal,drhox,dTx);
                mDy[na][nb] = extNSF_Durst::calcSelfDiffFlux(rhoVal,TVal,drhoy,dTy);
            }
        }

        //Order1 bat dau tu 1 vi khi order=0, div(phi)=0 nen tich phan bang 0
        for (int order1 = 1; order1 <= mathVar::orderElem; order1++)
        {
            double sumB(0.0), rhoOrder(rho[element][order1]);
            std::vector<std::vector<double>> A(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));

            //Avoid dividing by zero
            if (fabs(rhoOrder)<systemVar::epsilon)
            {
                if (rhoOrder>0)
                    rhoOrder=systemVar::epsilon;
                else
                    rhoOrder=-systemVar::epsilon;
            }

            //EcOrder phai nhan voi theta2 de tranh tich phan co gia tri unphysical o vi tri co strong discontinuity
            EcOrder=(-0.5*(rhou[element][order1]*rhou[element][order1] + rhov[element][order1]*rhov[element][order1])/rhoOrder)*theta2Arr[element];

            for (int na = 0; na <= mathVar::nGauss2D; na++)
            {
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                {
                    //Reset sumB
                    sumB=0.0;

                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                    for (int order2 = 0; order2 <= mathVar::orderElem; order2++)
                    {
                        if (elemType==3)
                            sumB += mathVar::BPts_Tri[order2][nanb];
                        else if (elemType==4)
                            sumB += mathVar::BPts_Quad[order2][nanb];
                    }
                    dBx = math::Calc_dBxdBy(element, order1, na, nb, 1);
                    dBy = math::Calc_dBxdBy(element, order1, na, nb, 2);

                    //Calculate Gauss values
                    A[na][nb]=EcOrder*sumB*(dBx*mDx[na][nb] + dBy*mDy[na][nb]);
                }
            }

            Int = math::volumeInte(A, element);

            EnergyEqnVolIntTerm[order1] = EnergyEqnVolIntTerm[order1] + Int;
        }
    }

    double calcSelfDiffFlux(double rho, double T, double gradRho, double gradT)
    {
        extNSF_Durst::applyBlendingFactorToDm();
        return (-extNSF_Durst::realDm*(gradRho/rho + gradT/(2*T)));
    }

    double calcDiffVelocity(std::vector<double> &U, std::vector<double> &dU)
    {
        double du(0.0), dv(0.0), dE(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0), rhoEVal(0.0), TVal;
        double drho(0.0),
            drhou(0.0),
            drhov(0.0),
            drhoE(0.0),
            dT(0.0);
        double vD(0.0);

        double uVal(0.0), vVal(0.0), eVal(0.0);

        rhoVal = U[0];
        rhouVal = U[1];
        rhovVal = U[2];
        rhoEVal = U[3];

        TVal = math::CalcTFromConsvVar(rhoVal,rhouVal,rhovVal,rhoEVal);
        eVal = material::Cv*TVal;

        uVal = rhouVal / rhoVal;
        vVal = rhovVal / rhoVal;

        drho = dU[0];
        drhou = dU[1];
        drhov = dU[2];
        drhoE = dU[3];

        du = math::calcRhouvEDeriv(drhou, drho, rhouVal, rhoVal);
        dv = math::calcRhouvEDeriv(drhov, drho, rhovVal, rhoVal);
        dE = math::calcRhouvEDeriv(drhoE, drho, rhoEVal, rhoVal);
        dT = math::calcTDeriv(dE, du, dv, uVal, vVal);

        vD = extNSF_Durst::calcSelfDiffFlux(rhoVal,TVal,drho,dT)/rhoVal;

        return vD;
    }
}

