#include "BCSupportFncs.h"
#include <vector>
#include <math.h>
#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGProcLib.h"
#include <algorithm>
#include <iostream>
#include "./parallelFunctions/cellData.h"
#include "./parallelFunctions/generalParallelFuncs.h"
#include "./parallelFunctions/parallelVariables.h"
//Boundary conditions header
#include "zeroGradient.h"
#include "symmetry.h"
#include "fixedValue.h"
#include "bcVariables.h"

//Durst model
#include "./extNSFEqns/FranzDurst/DurstModel.h"
#include "./extNSFEqns/FranzDurst/boundaryConditions.h"

//Custom boundary conditions
//Key word: Custom_Boundary_Conditions

//1. zeroRhoGradUncorectP
#include "./customBCs/zeroRhoGradUncorectP/p_zeroRhoGradUncorrectP.h"
#include "./customBCs/zeroRhoGradUncorectP/rho_zeroRhoGradUncorrectP.h"

//2. zeroRhoGradUncorectP
#include "./customBCs/reflectRhoGrad/p_reflectRhoGrad.h"
#include "./customBCs/reflectRhoGrad/rho_reflectRhoGrad.h"

//3. interiorSide
#include "./customBCs/interiorSide/p_interiorSide.h"
#include "./customBCs/interiorSide/T_interiorSide.h"
#include "./customBCs/interiorSide/u_interiorSide.h"
#include "./customBCs/interiorSide/rho_interiorSide.h"

//4. MaxwellSlip
#include "./customBCs/nonEquilibriumBCs/MaxwellSlip/u_MaxwellSlip.h"

//5. SmoluchowskyTJump
#include "./customBCs/nonEquilibriumBCs/SmoluchowskyTJump/T_SmoluchowskyTJump.h"

//6. zeroRhoGradUncorectP
#include "./customBCs/zeroRhoGrad/p_zeroRhoGrad.h"
#include "./customBCs/zeroRhoGrad/rho_zeroRhoGrad.h"

//7. Time varying BCs
#include "./boundaryConditions/customBCs/timeVaryingBCs/timeVaryingBCs_GenFuncs.h"
//7.1 waveTransmissive
#include "./boundaryConditions/customBCs/timeVaryingBCs/waveTransmissive/scalar_waveTransmissive.h"
#include "./boundaryConditions/customBCs/timeVaryingBCs/waveTransmissive/vector_waveTransmissive.h"

namespace BCSupportFncs
{
    void calcBCValuesAtUnfixedValueBCs()
    {
        int globleEdge(0);
        for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
        {
            globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
            int edgeGrp(auxUlti::getGrpOfEdge(globleEdge));
            int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);

            int element(0);
            std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(globleEdge);

            //Field U
            if (
                    //Cac dieu kien bien can tinh gia tri SurfaceBCFields tu conservative variables
                    UType==BCVars::velocityBCId::waveTransmissive //wave transmissive
                    )
            {
                for (int nG=0; nG<=mathVar::nGauss1D; nG++)
                {
                    double a,b;
                    std::tie(a, b) = auxUlti::getGaussSurfCoor(globleEdge, element, nG);
                    double rhoBC(math::pointValue(element,a,b,1,2)),
                            rhouBC(math::pointValue(element,a,b,2,2)),
                            rhovBC(math::pointValue(element,a,b,3,2));

                    SurfaceBCFields::uBc[ilocalEdge][nG]=rhouBC/rhoBC;
                    SurfaceBCFields::vBc[ilocalEdge][nG]=rhovBC/rhoBC;
                }
            }

            //Field T
            if (
                    //Cac dieu kien bien can tinh gia tri SurfaceBCFields tu conservative variables
                    TType==BCVars::temperatureBCId::waveTransmissive //wave transmissive
                    )
            {
                for (int nG=0; nG<=mathVar::nGauss1D; nG++)
                {
                    double a,b;
                    std::tie(a, b) = auxUlti::getGaussSurfCoor(globleEdge, element, nG);
                    double rhoBC(math::pointValue(element,a,b,1,2)),
                            rhouBC(math::pointValue(element,a,b,2,2)),
                            rhovBC(math::pointValue(element,a,b,3,2)),
                            rhoEBC(math::pointValue(element,a,b,4,2));

                    SurfaceBCFields::TBc[ilocalEdge][nG]=math::CalcTFromConsvVar(rhoBC,rhouBC,rhovBC,rhoEBC);
                }
            }

            //Field P
            if (
                    //Cac dieu kien bien can tinh gia tri SurfaceBCFields tu conservative variables
                    pType==BCVars::pressureBCId::waveTransmissive //wave transmissive
                    )
            {
                for (int nG=0; nG<=mathVar::nGauss1D; nG++)
                {
                    double a,b;
                    std::tie(a, b) = auxUlti::getGaussSurfCoor(globleEdge, element, nG);
                    double rhoBC(math::pointValue(element,a,b,1,2)),
                            rhouBC(math::pointValue(element,a,b,2,2)),
                            rhovBC(math::pointValue(element,a,b,3,2)),
                            rhoEBC(math::pointValue(element,a,b,4,2));

                    SurfaceBCFields::pBc[ilocalEdge][nG]=math::CalcP(math::CalcTFromConsvVar(rhoBC,rhouBC,rhovBC,rhoEBC),rhoBC);
                }
            }
        }
    }

    bool checkInflow(double u, double v, double nx, double ny)
    {
        bool inflow(true);
        double normUMag(0.0);
        std::vector < double > U(2, 0.0);
        std::vector < double > normVector(2, 0.0);

        U[0] = u;
        U[1] = v;
        normVector[0] = nx;
        normVector[1] = ny;

        normUMag = math::vectorDotProduct(U, normVector);

        if ((normUMag >= 0) || fabs(normUMag) <=0.00001)
        {
            inflow = false;
        }

        return inflow;
    }

    void decompseU(std::vector<double> &priVars, const std::vector<double> &U, bool calcT)
    {
        /* Ham decompose primary vars tu conservative vars (U). Trong truong ho
         * mass diffusion, U
         * Neu calcT = true, ham se giai T tu U va dU, sau do giai p tu T
         * Neu calcT = false, T va p = 0
        */

        double u(U[1]/U[0]),
                v(U[2]/U[0]),
                rho(U[0]),
                p(0.0), T(0.0);

        if (calcT)
        {
            T=math::CalcTFromConsvVar(U[0],U[1],U[2],U[3]);
            p=math::CalcP(T,rho);
        }

        priVars[0]=rho;
        priVars[1]=u;
        priVars[2]=v;
        priVars[3]=p;
        priVars[4]=T;
    }

    void decompsedU(std::vector<double> &priVars, const std::vector<double> &U, const std::vector<double> &dU, double T)
    {
        /* Ham decompose dao ham cua primary vars tu dao ham cua conservative vars (U). Trong truong co
         * mass diffusion, drhoE la dao ham cua convective energy.
        */

        double du, dv, dE, dT,
                drho(dU[0]), drhou(dU[1]), drhov(dU[2]), drhoE(dU[3]),
                rho(U[0]),rhou(U[1]),rhov(U[2]),rhoE(U[3]),
                mu(math::CalcVisCoef(T));

        //Remove mu
        drho=drho/mu;
        drhou=drhou/mu;
        drhov=drhov/mu;
        drhoE=drhoE/mu;

        du = math::calcRhouvEDeriv(drhou, drho, rhou, rho);
        dv = math::calcRhouvEDeriv(drhov, drho, rhov, rho);
        dE = math::calcRhouvEDeriv(drhoE, drho, rhoE, rho);
        dT = math::calcTDeriv(dE, du, dv, rhou/rho, rhov/rho);

        priVars[0]=drho;
        priVars[1]=du;
        priVars[2]=dv;
        priVars[3]=material::R*(dT*rho+drho*T);
        priVars[4]=dT;
    }

    void reconstructdU(std::vector<double> &dU, const std::vector<double> &priVars, const std::vector<double> &U, double T)
    {
        double drho(priVars[0]),
                du(priVars[1]),
                dv(priVars[2]),
                dT(priVars[4]), dE(0.0),
                rho(U[0]), u(U[1]/U[0]), v(U[2]/U[0]), E(U[3]/U[0]),
                mu(math::CalcVisCoef(T));

        dU[0]=drho*mu;
        dU[1]=(du*rho+drho*u)*mu;
        dU[2]=(dv*rho+drho*v)*mu;
        dE=dT*(material::R/(material::gamma-1))+u*du+v*dv;
        dU[3]=(dE*rho+drho*E)*mu;
    }

    void reconstructConvectiveU(std::vector<double> &U, const std::vector<double> &priVars)
    {
        double rho(priVars[0]),
                u(priVars[1]),
                v(priVars[2]),
                T(priVars[4]);
        U[0]=rho;
        U[1]=rho*u;
        U[2]=rho*v;
        U[3]=rho*(material::Cv*T+0.5*(u*u+v*v));
    }

    void correctPriVars(int edge, int edgeGrp, int nG, std::vector<double> &priVarsM, const std::vector<double> &priVarsP, const std::vector<double> &priVarsMean, const std::vector<double> &n, bool inflow)
    {
        /* priVars[0]=rho;
        priVars[1]=u;
        priVars[2]=v;
        priVars[3]=p;
        priVars[4]=T;
        */

        double //rhoP(priVarsP[0]),
                uP(priVarsP[1]),
                vP(priVarsP[2]),
                pP(priVarsP[3]),
                TP(priVarsP[4]),
                uMean(priVarsMean[1]),
                vMean(priVarsMean[2]),
                pMean(priVarsMean[3]),
                TMean(priVarsMean[4]);
        double rhoM(0.0), uM(0.0), vM(0.0), pM(0.0), TM(0.0);

        //Correct u, v, T, p
        //NOTE (22/08/2021): correct T truoc u de phu hop voi khi chay dieu kien bien nonequilibrium
        pM = correctPressure(edge,edgeGrp,nG,pP,pMean,inflow);
        TM = correctTemperature(edge,edgeGrp,nG,TP,TMean,inflow);
        std::tie(uM, vM) = correctVelocity(edge,edgeGrp,nG,uP,uMean,vP,vMean,n,inflow);

        priVarsM[0]=rhoM;
        priVarsM[1]=uM;
        priVarsM[2]=vM;
        priVarsM[3]=pM;
        priVarsM[4]=TM;

        //Vi ham correct rho lay gia tri tu priVarsM nen phai luu gia tri vao priVarsM roi moi chay ham correctRho
        BCSupportFncs::correctDensity(priVarsM,edgeGrp,priVarsP);

        //Bound T and p
        if (priVarsM[4]>limitVal::TUp)
        {
            priVarsM[4]=limitVal::TUp;
            //priVarsM[3]=math::CalcP(limitVal::TUp,priVarsM[0]);
        }
        else if (priVarsM[4]<limitVal::TDwn)
        {
            priVarsM[4]=limitVal::TDwn;
            //priVarsM[3]=math::CalcP(limitVal::TDwn,priVarsM[0]);
        }
    }

    void correctPriVarsGrad(int edge, int edgeGrp, int nG, std::vector<double> &dpriVarsXM, std::vector<double> &dpriVarsYM, const std::vector<double> &dpriVarsXP, const std::vector<double> &dpriVarsYP, const std::vector<double> &UP, const std::vector<double> &UM, double TP, double TM, const std::vector<double> &n, bool inflow)
    {
        /* priVars[0]=rho;
        priVars[1]=u;
        priVars[2]=v;
        priVars[3]=p;
        priVars[4]=T;
        */

        double drhoXP(dpriVarsXP[0]),
                duXP(dpriVarsXP[1]),
                dvXP(dpriVarsXP[2]),
                dpXP(dpriVarsXP[3]),
                dTXP(dpriVarsXP[4]),

                drhoYP(dpriVarsYP[0]),
                duYP(dpriVarsYP[1]),
                dvYP(dpriVarsYP[2]),
                dpYP(dpriVarsYP[3]),
                dTYP(dpriVarsYP[4]);

                /*
                drhoXMean(dpriVarsXMean[0]),
                duXMean(dpriVarsXMean[1]),
                dvXMean(dpriVarsXMean[2]),
                dpXMean(dpriVarsXMean[3]),
                dTXMean(dpriVarsXMean[4]),

                drhoYMean(dpriVarsYMean[0]),
                duYMean(dpriVarsYMean[1]),
                dvYMean(dpriVarsYMean[2]),
                dpYMean(dpriVarsYMean[3]),
                dTYMean(dpriVarsYMean[4]);
                        */
        double //rhoM(UM[0]),

                drhoXM(0.0),
                duXM(0.0),
                dvXM(0.0),
                dpXM(0.0),
                dTXM(0.0),

                drhoYM(0.0),
                duYM(0.0),
                dvYM(0.0),
                dpYM(0.0),
                dTYM(0.0);

        //NOTE (22/08/2021): correct T truoc u de phu hop voi khi chay dieu kien bien non equilibrium
        //Correct grad(T)
        std::tie(dTXM,dTYM)=correctTemperatureGrad(edgeGrp,nG,dTXP,dTYP,inflow,n);

        //Correct grad(u)
        std::tie(duXM,duYM)=correctVelocityGrad(edgeGrp,nG,duXP,duYP,inflow,n);
        std::tie(dvXM,dvYM)=correctVelocityGrad(edgeGrp,nG,dvXP,dvYP,inflow,n);

        //Correct grad(p)
        //Co the khong can thiet!
        std::tie(dpXM,dpYM)=correctPressureGrad(edgeGrp,dpXP,dpYP,inflow,n);

        //Correct grad(rho)
        std::tie(drhoXM, drhoYM)=correctDensityGrad(edgeGrp, dpXM, dpYM, dTXM, dTYM, drhoXP, drhoYP, UM, TM, n);

        //Save to output array
        dpriVarsXM[0]=drhoXM;
        dpriVarsXM[1]=duXM;
        dpriVarsXM[2]=dvXM;
        dpriVarsXM[3]=dpXM;
        dpriVarsXM[4]=dTXM;

        dpriVarsYM[0]=drhoYM;
        dpriVarsYM[1]=duYM;
        dpriVarsYM[2]=dvYM;
        dpriVarsYM[3]=dpYM;
        dpriVarsYM[4]=dTYM;
    }

    std::tuple <double, double> correctVelocity(int edge, int edgeGrp, int nG, double uP, double uMean, double vP, double vMean, const std::vector<double> &n, bool inflow)
    {
        //Ham correct gia tri bien U tai BC theo dieu kien bien
        /* Velocity boundary conditions:
            Id = 1:---------------------------------------------------------------------------------------------------------
                fixedValue
                value u v w
                Description: ham nay set dieu kien fixedValue cho patch.

            Id = 2:---------------------------------------------------------------------------------------------------------
                noSlip

            Id = 3:---------------------------------------------------------------------------------------------------------
                movingWall
                value u v w
                Description: set velocity cho wall (typically: cavity case)

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                inletValue u v w
                Description: Neu subsonic fixedValue neu inFlow va zeroGradient neu outFlow. Neu supersonic deu la zeroGradient.
                Nen dung cho outlet.

            Id = 5:---------------------------------------------------------------------------------------------------------
                slip
                sigmaU value
                uWall u v w
                Description: apply dieu kien slip cua Maxwell

            Id = 6:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry

            Id = 10:--------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song
        */

        int UType(bcValues::UBcType[edgeGrp - 1]), //Lay UType de xac dinh kieu dk bien
            loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::vector<double> UP{uP,vP},UMean{uMean,vMean},
        UM(2, 0.0),
        UBC{SurfaceBCFields::uBc[loc][nG],SurfaceBCFields::vBc[loc][nG]},
        UZero{0.0,0.0};

        bool isStrong(BCVars::DirichletAppMethUStrong[edgeGrp -1]);

        switch (UType)
        {
        case BCVars::velocityBCId::fixedValue:
            {
                fixedValue_vector(UM,UP,UBC,isStrong);
            }
            break;
        case BCVars::velocityBCId::noSlip:
            {
                fixedValue_vector(UM,UP,UZero,isStrong);
            }
            break;
        case BCVars::velocityBCId::movingWall:
            {
                fixedValue_vector(UM,UP,UBC,isStrong);
            }
            break;
        case BCVars::velocityBCId::inletOutlet:
            {
                if (flowProperties::subsonic)
                {
                    if (inflow)
                    {
                        fixedValue_vector(UM,UP,UBC,isStrong);
                    }
                    else
                    {
                        zeroGradient_vector(UM,UP);
                    }
                }
                else
                {
                    zeroGradient_vector(UM,UP);
                }
            }
            break;
        case BCVars::velocityBCId::slipWall:
            {
                fixedValue_vector(UM,UP,UBC,isStrong);
            }
            break;
        case BCVars::velocityBCId::zeroGrad:
            {
                zeroGradient_vector(UM,UP);
            }
            break;

//------//Custom_Boundary_Conditions----------------------------------------
        //interiorSide
        case BCVars::velocityBCId::interiorSide:
            {
                interiorSide::correctU(UM,UP);
            }
            break;
        //Maxwell Slip
        case BCVars::velocityBCId::MaxwellSlip:
            {
                MaxwellSlip::correctU(edge,edgeGrp,nG,UM,UP);
            }
            break;
        //Wave transmissive - 07/06/2022
        case BCVars::velocityBCId::waveTransmissive:
            {
                //Tuong tu fixed value
                fixedValue_vector(UM,UP,UBC,isStrong);
            }
            break;
//------//Custom_Boundary_Conditions----------------------------------------

        case BCVars::generalBCId::matched:
            {
                //matched BC
            }
            break;
            default:
            {
                std::cout<<"Warning! Unknow velocity boudary type Id = "<<UType<<" .\n";
            }
            break;
        }

        return std::make_tuple(UM[0],UM[1]);
    }

    std::tuple <double, double> correctVelocityGrad(int edgeGrp, int nG, double duXP, double duYP, bool inflow, const std::vector<double> &n)
    {
        //Ham correct gia tri bien U tai BC theo dieu kien bien
        /* Velocity boundary conditions:
            Id = 1:---------------------------------------------------------------------------------------------------------
                fixedValue
                value u v w
                Description: ham nay set dieu kien fixedValue cho patch.

            Id = 2:---------------------------------------------------------------------------------------------------------
                noSlip

            Id = 3:---------------------------------------------------------------------------------------------------------
                movingWall
                value u v w
                Description: set velocity cho wall (typically: cavity case)

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                inletValue u v w
                Description: Neu subsonic fixedValue neu inFlow va zeroGradient neu outFlow. Neu supersonic deu la zeroGradient.
                Nen dung cho outlet.

            Id = 5:---------------------------------------------------------------------------------------------------------
                slip
                sigmaU value
                uWall u v w
                Description: apply dieu kien slip cua Maxwell

            Id = 6:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry

            Id = 10:--------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song
        */

        int UType(bcValues::UBcType[edgeGrp - 1]); //Lay UType de xac dinh kieu dk bien
        std::vector<double> duP{duXP, duYP}, duM(2, 0.0);

        bool isStrong(BCVars::NewmannAppMethGradUStrong[edgeGrp -1]);

        switch (UType)
        {
            case BCVars::velocityBCId::fixedValue:
            {
                duM=duP;
            }
            break;
        case BCVars::velocityBCId::noSlip:
            {
                duM=duP;
            }
            break;
        case BCVars::velocityBCId::movingWall:
            {
                duM=duP;
            }
            break;
        case BCVars::velocityBCId::inletOutlet:
            {
                if (flowProperties::subsonic)
                {
                    if (inflow)
                    {
                        duM=duP;
                    }
                    else
                    {
                        zeroGradient_correctGrad(duM,duP,n,isStrong);
                    }
                }
                else
                {
                    zeroGradient_correctGrad(duM,duP,n,isStrong);
                }
            }
            break;
        case BCVars::velocityBCId::slipWall:
            {
                duM=duP;
            }
            break;
        case BCVars::velocityBCId::zeroGrad:
            {
                zeroGradient_correctGrad(duM,duP,n,isStrong);
            }
            break;

//------//Custom_Boundary_Conditions----------------------------------------
        //interiorSide
        case BCVars::velocityBCId::interiorSide:
            {
                interiorSide::correctGradU(duM,duP);
            }
            break;
        //Maxwell Slip
        case BCVars::velocityBCId::MaxwellSlip:
            {
                //Giong fixedValue
                duM=duP;
            }
            break;
        //Wave transmissive - 07/06/2022
        case BCVars::velocityBCId::waveTransmissive:
            {
                //Tuong tu fixed value
                duM=duP;
            }
            break;
//------//Custom_Boundary_Conditions----------------------------------------

        case BCVars::generalBCId::matched:
            {
                //matched BC
            }
            break;
            default:
            {
                std::cout<<"Warning! Unknow velocity boudary type Id = "<<UType<<" .\n";
            }
            break;
        }

        return std::make_tuple(duM[0],duM[1]);
    }

    double correctTemperature(int edge, int edgeGrp, int nG, double TP, double TMean, bool inflow)
    {
        //Ham correct gia tri bien T tai BC theo dieu kien bien
        /* Temperature boundary conditions:
            Id = 2:---------------------------------------------------------------------------------------------------------
                fixedValue
                value T
                ==> grad surface = grad surface plus

            Id = 3:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid. Argument TP de nguyen, de phong truong hop
                muon lay zeroGradient la gia tri tai diem Gauss dang tinh cua cell phia plus.
                ==> giai zero normal grad

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                value T
                Description:
                ==> giai zero normal grad neu outflow, grad surface = grad mean neu inflow

            Id = 6:---------------------------------------------------------------------------------------------------------
                temperatureJump
                sigmaT value
                TWall T
                Description: apply dieu kien temperature cua Smoluchowsky
                ==> giai zero normal grad

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry
                ==> reflect vector grad T

            Id = 10:---------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song
        */

        int TType(bcValues::TBcType[edgeGrp - 1]), //Lay UType de xac dinh kieu dk bien
            loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        double TBC(SurfaceBCFields::TBc[loc][nG]), TM(0);

        bool isStrong(BCVars::DirichletAppMethTStrong[edgeGrp-1]);

        switch (TType)
        {
        case BCVars::temperatureBCId::fixedValue:
            {
                TM = fixedValue_scalar(TP,TBC,isStrong);
            }
            break;
        case BCVars::temperatureBCId::zeroGrad:
            {
                TM = zeroGradient_scalar(TP);
            }
            break;
        case BCVars::temperatureBCId::inletOutlet:
            if (flowProperties::subsonic)
            {
                if (inflow)
                {
                    TM = fixedValue_scalar(TP,TBC,isStrong);
                }
                else
                {
                    TM = zeroGradient_scalar(TP);
                }
            }
            else
            {
                TM = zeroGradient_scalar(TP);
            }
            break;
        case BCVars::temperatureBCId::temperatureJump:
            {
                TM = fixedValue_scalar(TP,TBC,isStrong);
            }
            break;
//------//Custom_Boundary_Conditions----------------------------------------
        //interiorSide
        case BCVars::temperatureBCId::interiorSide:
            {
                interiorSide::correctT(TM,TP);
            }
            break;
        //Smoluchowsky Jump
        case BCVars::temperatureBCId::SmoluchowskyTJump:
            {
                SmoluchowskyTJump::correctT(edge,edgeGrp,nG,TM,TP);
            }
            break;
        //Wave transmissive - 07/06/2022
        case BCVars::temperatureBCId::waveTransmissive:
            {
                //Tuong tu fixed value
                TM = fixedValue_scalar(TP,TBC,isStrong);
            }
            break;
//------//Custom_Boundary_Conditions----------------------------------------
        case BCVars::generalBCId::matched:
            {
                //matched BC
            }
            break;
        default:
            {
                std::cout<<"Warning! Unknow temperature boudary type Id = "<<TType<<" .\n";
            }
            break;
        }

        return TM;
    }

    std::tuple<double, double> correctTemperatureGrad(int edgeGrp, int nG, double dTXP, double dTYP, bool inflow, const std::vector<double> &n)
    {
        //Ham correct gia tri bien dT tai BC theo dieu kien bien.
        /* Temperature boundary conditions:
            Id = 2:---------------------------------------------------------------------------------------------------------
                fixedValue
                value T
                ==> grad surface = grad surface plus

            Id = 3:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid. Argument TP de nguyen, de phong truong hop
                muon lay zeroGradient la gia tri tai diem Gauss dang tinh cua cell phia plus.
                ==> giai zero normal grad

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                value T
                Description:
                ==> giai zero normal grad neu outflow, grad surface = grad mean neu inflow

            Id = 6:---------------------------------------------------------------------------------------------------------
                temperatureJump
                sigmaT value
                TWall T
                Description: apply dieu kien temperature cua Smoluchowsky
                ==> grad surface = grad surface plus

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry
                ==> reflect vector grad T

            Id = 10:---------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song
        */

        int TType(bcValues::TBcType[edgeGrp - 1]); //Lay UType de xac dinh kieu dk bien
        std::vector<double> dTP{dTXP, dTYP}, dTM(2, 0.0);

        bool isStrong(BCVars::NewmannAppMethGradTStrong[edgeGrp-1]);

        switch (TType)
        {
        case BCVars::temperatureBCId::fixedValue:
            {
                dTM=dTP;
            }
            break;
        case BCVars::temperatureBCId::zeroGrad:
            {
                zeroGradient_correctGrad(dTM,dTP,n,isStrong);
            }
            break;
        case BCVars::temperatureBCId::inletOutlet:
            if (flowProperties::subsonic)
            {
                if (inflow)
                {
                    dTM=dTP;
                }
                else
                {
                    zeroGradient_correctGrad(dTM,dTP,n,isStrong);
                }
            }
            else
            {
                zeroGradient_correctGrad(dTM,dTP,n,isStrong);
            }
            break;
        case BCVars::temperatureBCId::temperatureJump:
            {
                dTM=dTP;
            }
            break;

//------//Custom_Boundary_Conditions----------------------------------------
        //interiorSide
        case BCVars::temperatureBCId::interiorSide:
            {
                interiorSide::correctGradT(dTM,dTP);
            }
            break;
        //Smoluchowsky Jump
        case BCVars::temperatureBCId::SmoluchowskyTJump:
            {
                //Giong fixedValue
                dTM=dTP;
            }
            break;
        //Wave transmissive - 07/06/2022
        case BCVars::temperatureBCId::waveTransmissive:
            {
                //Tuong tu fixed value
                dTM=dTP;
            }
            break;
//------//Custom_Boundary_Conditions----------------------------------------

        case BCVars::generalBCId::matched:
            {
                //matched BC
            }
            break;
        default:
            {
                std::cout<<"Warning! Unknow temperature boudary type Id = "<<TType<<" .\n";
            }
            break;
        }

        return std::make_tuple(dTM[0],dTM[1]);
    }

    double correctPressure(int edge, int edgeGrp, int nG, double pP, double pMean, bool inflow)
    {
        //Ham correct gia tri bien P tai BC theo dieu kien bien
        /* Pressure boundary conditions:
            Id = 1:---------------------------------------------------------------------------------------------------------
                fixedValue
                value p

            Id = 2:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid. Argument pP de nguyen, de phong truong hop
                muon lay zeroGradient la gia tri tai diem Gauss dang tinh cua cell phia plus.

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                value p
                Description:

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry

            Id = 10:---------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song

            Id = 11:---------------------------------------------------------------------------------------------------------
                interpFrmDensity
                Description: Dung cho wall. Gia tri cua P duoc tinh tu rhoP va Twall, van dung zeroGradient cho P khi tinh gradient.
                correctRho = true khi dung kieu BC nay.
        */

        int pType(bcValues::pBcType[edgeGrp - 1]), //Lay UType de xac dinh kieu dk bien
            loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        double pBC(SurfaceBCFields::pBc[loc][nG]), pM(0);

        bool isStrong(BCVars::DirichletAppMethPStrong[edgeGrp-1]);

        switch (pType)
        {
        case BCVars::pressureBCId::fixedValue:
            {
                pM = fixedValue_scalar(pP,pBC,isStrong);
            }
            break;
        case BCVars::pressureBCId::zeroGrad:
            {
                pM = zeroGradient_scalar(pP);
            }
            break;
        case BCVars::pressureBCId::inletOutlet:
            {
                if (flowProperties::subsonic)
                {
                    if (inflow)
                    {
                        pM = fixedValue_scalar(pP,pBC,isStrong);
                    }
                    else
                    {
                        pM = zeroGradient_scalar(pP);
                    }
                }
                else
                {
                    pM = zeroGradient_scalar(pP);
                }
            }
            break;
        case BCVars::generalBCId::matched:
            {
                //matched BC
            }
            break;
        case BCVars::pressureBCId::interpFrmDensity:
            {
                pM = zeroGradient_scalar(pP); //==> update lai p theo TM va rhoP o ngoai ham nay
            }
            break;

//------//Custom_Boundary_Conditions----------------------------------------
        //zeroGradRhoUncorrectP
        case BCVars::pressureBCId::zeroGradRhoUncorrectP:
            {
                zeroRhoGradUncorectP::correctP(pM,pP);
            }
            break;
            //zeroGradRhoUncorrectP
        case BCVars::pressureBCId::reflectGradRho:
            {
                reflectRhoGrad::correctP(pM,pP);
            }
            break;
        case BCVars::pressureBCId::interiorSide:
            {
                interiorSide::correctP(pM,pP);
            }
            break;
        case BCVars::pressureBCId::zeroRhoGrad:
            {
                zeroRhoGrad::correctP(pM,pP);
            }
            break;

        //Wave transmissive - 07/06/2022
        case BCVars::pressureBCId::waveTransmissive:
            {
                //Tuong tu fixed value
                pM = fixedValue_scalar(pP,pBC,isStrong);
            }
            break;
//------//Custom_Boundary_Conditions----------------------------------------


        default:
            {
                std::cout<<"Warning! Unknow pressure boudary type Id = "<<pType<<".\n";
            }
            break;
        }

        return pM;
    }

    std::tuple<double, double> correctPressureGrad(int edgeGrp, double dpXP, double dpYP, bool inflow, const std::vector<double> &n)
    {
        //Ham correct gia tri bien P tai BC theo dieu kien bien
        /* Pressure boundary conditions:
            Id = 1:---------------------------------------------------------------------------------------------------------
                fixedValue
                value p
                ==> grad surface = grad surface plus

            Id = 2:---------------------------------------------------------------------------------------------------------
                zeroGradient
                Description: zeroGradient lay theo gia tri trung binh tai cell centroid. Argument pP de nguyen, de phong truong hop
                muon lay zeroGradient la gia tri tai diem Gauss dang tinh cua cell phia plus.
                ==> giai zero normal grad

            Id = 4:---------------------------------------------------------------------------------------------------------
                inletOutlet
                value p
                Description:
                ==> giai zero normal grad neu outflow, grad surface = grad mean neu inflow

            Id = 7:---------------------------------------------------------------------------------------------------------
                symmetry
                ==> reflect vector grad T

            Id = 10:---------------------------------------------------------------------------------------------------------
                matched
                Description: dung cho tinh toan song song

            Id = 11:---------------------------------------------------------------------------------------------------------
                interpFrmDensity
                Description: Dung cho wall. Gia tri cua P duoc tinh tu rhoP va Twall, van dung zeroGradient cho P khi tinh gradient.
                correctRho = true khi dung kieu BC nay.
        */

        int pType(bcValues::pBcType[edgeGrp - 1]); //Lay UType de xac dinh kieu dk bien
        std::vector<double> dpP{dpXP, dpYP}, dpM(2, 0.0);

        bool isStrong(BCVars::NewmannAppMethGradPStrong[edgeGrp-1]);

        switch (pType)
        {
        case BCVars::pressureBCId::fixedValue:
            {
                dpM=dpP;
            }
            break;
        case BCVars::pressureBCId::zeroGrad:
            {
                zeroGradient_correctGrad(dpM,dpP,n,isStrong);
            }
            break;
        case BCVars::pressureBCId::inletOutlet:
            {
                if (flowProperties::subsonic)
                {
                    if (inflow)
                    {
                        dpM=dpP;
                    }
                    else
                    {
                        zeroGradient_correctGrad(dpM,dpP,n,isStrong);
                    }
                }
                else
                {
                    zeroGradient_correctGrad(dpM,dpP,n,isStrong);
                }
            }
            break;
        case BCVars::generalBCId::matched:
            {
                //matched BC
            }
            break;
        case BCVars::pressureBCId::interpFrmDensity:
            {
                //Test BC
                zeroGradient_correctGrad(dpM,dpP,n,isStrong);
            }
            break;


//------//Custom_Boundary_Conditions----------------------------------------
        //zeroGradRhoUncorrectP
        case BCVars::pressureBCId::zeroGradRhoUncorrectP:
            {
                zeroRhoGradUncorectP::correctGradP(dpM,dpP);
            }
            break;
        //reflectGradRho
        case BCVars::pressureBCId::reflectGradRho:
            {
                reflectRhoGrad::correctGradP(dpM,dpP);
            }
            break;
        //interiorSide
        case BCVars::pressureBCId::interiorSide:
            {
                interiorSide::correctGradP(dpM,dpP);
            }
            break;
        //interiorSide
        case BCVars::pressureBCId::zeroRhoGrad:
            {
                zeroRhoGrad::correctGradP(dpM,dpP);
            }
            break;
        //Wave transmissive - 07/06/2022
        case BCVars::pressureBCId::waveTransmissive:
            {
                //Tuong tu fixed value
                dpM=dpP;
            }
            break;
//------//Custom_Boundary_Conditions----------------------------------------


        default:
            {
                std::cout<<"Warning! Unknow pressure boudary type Id = "<<pType<<" .\n";
            }
            break;
        }

        return std::make_tuple(dpM[0], dpM[1]);
    }

    void correctDensity(std::vector<double> &priVarsM, int edgeGrp, const std::vector<double> &priVarsP)
    {
        double rhoP(priVarsP[0]);

        //General: Correct rhoM theo pM va TM theo ideal gas
        priVarsM[0]=priVarsM[3]/(material::R*priVarsM[4]);

        //Correct rho theo tung dieu kien rieng biet
        int pType(bcValues::pBcType[edgeGrp - 1]);
        if (pType==BCVars::pressureBCId::interpFrmDensity)
        {
            priVarsM[0]=rhoP;
        }

//------//Custom_Boundary_Conditions----------------------------------------
        //zeroGradRhoUncorrectP
        else if (pType==BCVars::pressureBCId::zeroGradRhoUncorrectP)
        {
            zeroRhoGradUncorectP::correctRho(priVarsM[0],rhoP);
        }
        //interiorSide
        else if (pType==BCVars::pressureBCId::interiorSide)
        {
            interiorSide::correctRho(priVarsM[0],rhoP);
        }
        //zeroRhoGrad
        else if (pType==BCVars::pressureBCId::zeroRhoGrad)
        {
            zeroRhoGrad::correctRho(priVarsM[0],rhoP);
        }
//------//Custom_Boundary_Conditions----------------------------------------
    }

    std::tuple<double, double> correctDensityGrad(int edgeGrp, double dpMX, double dpMY, double dTMX, double dTMY, double dRhoPX, double dRhoPY, const std::vector<double> &UM, double TM, const std::vector<double> &n)
    {
        std::vector<double> drhoP{dRhoPX, dRhoPY}, drhoM(2, 0.0);

        //General, dRho- = dRho+ (don't follow perfect gas equation)
        drhoM[0] = dRhoPX;
        drhoM[1] = dRhoPY;
        /*
        drhoM[0]=((dpMX/material::R)-dTMX*UM[0])/TM;
        drhoM[1]=((dpMY/material::R)-dTMY*UM[0])/TM;*/

        //Modify grad(rho) theo tung dieu kien rieng biet
        int pType(bcValues::pBcType[edgeGrp - 1]);
        bool isStrong(BCVars::NewmannAppMethGradPStrong[edgeGrp-1]);

//------//Custom_Boundary_Conditions----------------------------------------
        //zeroGradRhoUncorrectP
        if (pType==BCVars::pressureBCId::zeroGradRhoUncorrectP)
        {
            zeroRhoGradUncorectP::correctGradRho(drhoM, drhoP);
        }
        //reflectGradRho
        else if (pType==BCVars::pressureBCId::reflectGradRho)
        {
            reflectRhoGrad::correctGradRho(edgeGrp, drhoM, drhoP, n);
        }
        //interiorSide
        else if (pType==BCVars::pressureBCId::interiorSide)
        {
            interiorSide::correctGradRho(drhoM, drhoP);
        }
        //zeroRhoGrad
        else if (pType==BCVars::pressureBCId::zeroRhoGrad)
        {
            zeroRhoGrad::correctGradRho(drhoM, drhoP, n, isStrong);
        }
//------//Custom_Boundary_Conditions----------------------------------------

        return std::make_tuple(drhoM[0], drhoM[1]);
    }

    namespace auxilaryBCs {
        void getUPlus(int edge, int nG, std::vector<double> &UPlus)
        {
            UPlus[0] = surfaceFields::rho[edge][nG];
            UPlus[1] = surfaceFields::rhou[edge][nG];
            UPlus[2] = surfaceFields::rhov[edge][nG];
            UPlus[3] = surfaceFields::rhoE[edge][nG];
        }

        /*Ham apply dieu kien back flow tai vi tri bi back (reversed) flow*/
        void calcUReversedFlow(int edge, int nG, std::vector<double> &U)
        {
            int edgeGrp(auxUlti::getGrpOfEdge(edge));
            U[0] = (bcValues::pBCFixed[edgeGrp - 1] / (material::R*bcValues::TBCFixed[edgeGrp - 1]));
            U[1] = U[0] * bcValues::uBCFixed[edgeGrp - 1];
            U[2] = U[0] * bcValues::vBCFixed[edgeGrp - 1];
            U[3] = U[0] * (bcValues::TBCFixed[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBCFixed[edgeGrp - 1], 2) + pow(bcValues::vBCFixed[edgeGrp - 1], 2)));
            surfaceFields::T[edge][nG]=bcValues::TBCFixed[edgeGrp - 1];
        }

        void calcUMean(int element, std::vector<double> &UMean)
        {
            double TMean(0.0);
            auxUlti::getVectorUMeanOfCell(element, UMean);

            TMean = math::calcMeanPriVar(element,1);

            UMean[3]=UMean[0]*material::Cv*TMean+0.5*(UMean[1]*UMean[1]+UMean[2]*UMean[2])/UMean[0];
        }
    }

    namespace NSFEqBCs {
    void findMaxLxFConstantOnEdge(int element, int edge)
    {
        /* Ham tinh he so C tren 1 edge trong truong hop Flux type la LxF,
         * cell id input vao ham la master cell cua edge
        */

        double TM, TP;
        std::vector<double> UM(4), UP(4);
        std::vector<double> CArray(mathVar::nGauss1D+1);
        for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
        {
            TP=surfaceFields::T[edge][nG];
            TM=surfaceFields::T[edge][nG+mathVar::nGauss1D+1];

            for (int i = 0; i < 4; i++)
            {
                //std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
                std::tie(UP[i], UM[i]) = auxUlti::getUAtInterfaces(edge, element, nG, i + 1);
            }

            double uMagP(sqrt(pow(UP[1]/UP[0],2)+pow(UP[2]/UP[0],2))),
                    uMagM(sqrt(pow(UM[1]/UM[0],2)+pow(UM[2]/UM[0],2))),
                    aP(math::CalcSpeedOfSound(TP)),
                    aM(math::CalcSpeedOfSound(TM));

            CArray[nG]=(math::numericalFluxes::constantC(uMagP,uMagM,aP,aM));
        }
        LxFConst[edge]=*max_element(CArray.begin(), CArray.end());
        //-----------------------------------------------
    }

    void getdUPlus(int edge, int nG, std::vector<double> &dUXPlus, std::vector<double> &dUYPlus)
    {
        //int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));

        dUXPlus[0] = surfaceFields::dRhoX[edge][nG];
        dUXPlus[1] = surfaceFields::dRhouX[edge][nG];
        dUXPlus[2] = surfaceFields::dRhovX[edge][nG];
        dUXPlus[3] = surfaceFields::dRhoEX[edge][nG];

        dUYPlus[0] = surfaceFields::dRhoY[edge][nG];
        dUYPlus[1] = surfaceFields::dRhouY[edge][nG];
        dUYPlus[2] = surfaceFields::dRhovY[edge][nG];
        dUYPlus[3] = surfaceFields::dRhoEY[edge][nG];
    }

    void calcdUMean(int edge, int element, std::vector<double> &dUXMean, std::vector<double> &dUYMean)
    {
        //Compute dU+
        /*
        for (int i = 0; i < 4; i++)
        {
            dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
            dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
        }
        */
        dUXMean=math::pointSVars(edge,element,0,0,1,2);
        dUYMean=math::pointSVars(edge,element,0,0,2,2);
    }

    void NSFEqFluxes(int edgeId, std::vector<double> &Fluxes, double TPlus, double TMinus, std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector)
    {
        /*StressHeat matrix has form:
        [tauXx		tauXy		Qx]
        [tauYx		tauYy		Qy]
        */
        std::vector<std::vector<double>> StressHeatP(2, std::vector<double>(3, 0.0));
        std::vector<std::vector<double>> StressHeatM(2, std::vector<double>(3, 0.0));

        //Local coordinates. _M is mass diffusion
        std::vector<double> ULPlus(4, 0.0), ULMinus(4, 0.0);

        std::vector<double>
                //Local fluxes
                invisFluxLocalXPlus(4,0.0),
                invisFluxLocalYPlus(4,0.0),
                visFluxLocalXPlus(4,0.0),
                visFluxLocalYPlus(4,0.0),

                invisFluxLocalXMinus(4,0.0),
                invisFluxLocalYMinus(4,0.0),
                visFluxLocalXMinus(4,0.0),
                visFluxLocalYMinus(4,0.0),

                //vector lay gia tri flux o he toa do global
                invisFluxGlobal(4,0.0),
                visFluxGlobal(4,0.0);

        //double umP, vmP;

        //Correct rhoEM = rhoEP following Mengaldo
        //double rhoEM_org(UMinus[4]); //Backup
        //UMinus[4]=UPlus[4];

        /* Voi ham calcLocalInviscidFlux, neu flux type la LxF hoac central, cac bien
         * invisFluxLocalXMaster,.., ULMaster,.. deu luu gia tri tai he toa do global.
         * Neu flux type la Roe hoac HLL, cac bien nay luu gia tri tai he local.
        */

        process::NSFEq::calcLocalInviscidFlux(UPlus,invisFluxLocalXPlus,invisFluxLocalYPlus,ULPlus,normVector,TPlus);
        process::NSFEq::calcLocalInviscidFlux(UMinus,invisFluxLocalXMinus,invisFluxLocalYMinus,ULMinus,normVector,TMinus);


        //Recover rhoEM before calculating viscous term
        //UMinus[4]=rhoEM_org;

        /* Voi ham calcLocalViscousFlux, vi flux type luon la central, cac bien
         * invisFluxLocalXMaster,.., ULMaster,.. deu luu gia tri tai he toa do global.
        */
        if (flowProperties::viscous)
        {
            process::NSFEq::calcLocalViscousFlux(UPlus,dUXPlus,dUYPlus,visFluxLocalXPlus,visFluxLocalYPlus,normVector);
            process::NSFEq::calcLocalViscousFlux(UMinus,dUXMinus,dUYMinus,visFluxLocalXMinus,visFluxLocalYMinus,normVector);
        }

        //Tinh flux
        //Vi dang tinh cho cell master cua edge nen phia trong (phia +) la cua cell master
        math::numericalFluxes::calConvectiveFluxes(edgeId, invisFluxGlobal,
                                                   invisFluxLocalXPlus, invisFluxLocalXMinus,
                                                   invisFluxLocalYPlus, invisFluxLocalYMinus,
                                                   ULPlus, ULMinus,
                                                   TPlus, TMinus,
                                                   normVector);

        //Viscous flux tinh bang central flux
        if (flowProperties::viscous)
        {
            math::numericalFluxes::centralFlux(visFluxGlobal,
                                               visFluxLocalXPlus, visFluxLocalXMinus,
                                               visFluxLocalYPlus, visFluxLocalYMinus,
                                               normVector);
        }

        //Luu flux vao vector
        for (int i=0; i<4; i++)
        {
            Fluxes[i]=invisFluxGlobal[i]+visFluxGlobal[i];
        }

        /* Save mass diffusion velocity at wall
         * Use central flux
        */
        int localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edgeId));
        extNSF_Durst::uMassWall[localEdgeId]=0.5*(visFluxLocalXPlus[0] + visFluxLocalXMinus[0])/(0.5*(ULPlus[0]+ULMinus[0]));
        extNSF_Durst::vMassWall[localEdgeId]=0.5*(visFluxLocalYPlus[0] + visFluxLocalYMinus[0])/(0.5*(ULPlus[0]+ULMinus[0]));
    }
    }

    namespace parallel {

    void getUMinus_convective(int edge, int nG, std::vector<double> &UMinus)
    {
        //int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));

        UMinus[0] = surfaceFields::rho[edge][mathVar::nGauss1D+nG+1];
        UMinus[1] = surfaceFields::rhou[edge][mathVar::nGauss1D+nG+1];
        UMinus[2] = surfaceFields::rhov[edge][mathVar::nGauss1D+nG+1];
        UMinus[3] = surfaceFields::rhoE[edge][mathVar::nGauss1D+nG+1];
    }

    void getdUMinus(int edge, int nG, std::vector<double> &dUMinus, int dir)
    {
        //int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));

        if (dir==1)
        {
            dUMinus[0] = surfaceFields::dRhoX[edge][mathVar::nGauss1D+nG+1];
            dUMinus[1] = surfaceFields::dRhouX[edge][mathVar::nGauss1D+nG+1];
            dUMinus[2] = surfaceFields::dRhovX[edge][mathVar::nGauss1D+nG+1];
            dUMinus[3] = surfaceFields::dRhoEX[edge][mathVar::nGauss1D+nG+1];
        }
        else
        {
            dUMinus[0] = surfaceFields::dRhoY[edge][mathVar::nGauss1D+nG+1];
            dUMinus[1] = surfaceFields::dRhouY[edge][mathVar::nGauss1D+nG+1];
            dUMinus[2] = surfaceFields::dRhovY[edge][mathVar::nGauss1D+nG+1];
            dUMinus[3] = surfaceFields::dRhoEY[edge][mathVar::nGauss1D+nG+1];
        }
    }
    }
}
