#include "massDiffusion.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include <algorithm>
#include "DGPostProcessLib.h"
#include <math.h>
#include "DGProcLib.h"
#include "./boundaryConditions/DGBCsLib.h"
#include <vector>
#include <iostream>
#include "./limiters/detectTroubleCell.h"
#include "./limiters/limiterController.h"
#include "./limiters/parallelFuncs.h"
#include "./limiters/positivityPreserving/positivityPreserving.h"

#include "DGIOLib.h"

//Parallel functions
#include "./parallelFunctions/generalParallelFuncs.h"
#include "./parallelFunctions/GaussPointData.h"

#include "debuggingFuncs.h"

namespace limiter
{
    namespace massDiffusion
    {
        //Variables -----------------------------------------------------------------
        double DmCoeff(1.32),
            pseudoCo(0.001);
        int maxIter(50);

        bool trbCellAtMatchedBC(false);
        bool *markerOfTrbCellAtMatchedBC = new bool[1];
        bool *markerOfTrbCellAtMatchedBC_buffer = new bool[1];

        int numOfTrbEdge(0);

        //Debug
        int RKOrder(0); //DEBUG_COMM


        //Functions -----------------------------------------------------------------
        void limiter()
        {
            bool limit(false);
            int iter(0), totalNumTrbCells(0);

            //save current time step
            double dt_temp(dt);

            //modify dt
            dt*=limiter::massDiffusion::pseudoCo;

            limit = limiter::massDiffusion::checkTroubleCells();

            if (limit)
            {
                bool runFlag(true);
                do {
                    //Synch cell status
                    //limiter_parallelFuncs::massDiff::synchCellStatus();

                    //Giai pt NSF mo rong bang TVDRK3
                    for (int iRKOrder = 1; iRKOrder <= 3; iRKOrder++)
                    {
                        limiter::massDiffusion::RKOrder=iRKOrder; //DEBUG_COMM

                        limiter::massDiffusion::mathFuncs::limiter_1step(iRKOrder);
                        limiter::massDiffusion::updateVariables();
                        //limiter::massDiffusion::limitRho_PPLimitiedVer();
                    }

                    //Reset markerOfTrbCellAtMatchedBC
                    //auxUlti::initialize1DBoolArray(limiter::massDiffusion::markerOfTrbCellAtMatchedBC, meshVar::numBCEdges, false);

                    std::tie(runFlag, totalNumTrbCells)=limiter::massDiffusion::checkRunning();
                    iter++;
                }
                while (runFlag&&iter<=limiter::massDiffusion::maxIter);

                int maxIter(parallelFuncs_Gen::maxIntOverProcs(iter));
                if (systemVar::currentProc==0)
                {
                    std::cout <<"Final no of cells = "<<totalNumTrbCells<<", No of Iter: "<<maxIter<<".\n";
                }
            }

            dt=dt_temp;
        }

        /**
         * @brief Function detects trouble cells using PerssonPeraireDetector
         *
         * and marks on limitVal::troubleCellsMarker.
         *
         * @return
         */
        bool checkTroubleCells()
        {
            bool trouble(false);
            int counter(0), totalTrbCell(0);
            for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
            {
                if (troubledCellDetector::PerssonPeraireDetector(nelem))
                {
                    limitVal::troubleCellsMarker[nelem]=true;
                    counter++;
                    trouble=true;

                    //Check xem co trouble cell co edge la bien match hay khong, neu co thi moi active cac ham synch data
                    /*
                    int elemType(auxUlti::checkType(nelem)), edgeId(-1);
                    for (int i=0; i<elemType; i++)
                    {
                        edgeId = meshVar::inedel[nelem][i];
                        if (auxUlti::getBCType(i)==4) //Neu type = 4 thi la bien matched
                        {
                            //Mark edge cua trouble cell
                            limiter::massDiffusion::markerOfTrbCellAtMatchedBC[auxUlti::getAdressOfBCEdgesOnBCValsArray(edgeId)]=true;
                        }
                        else
                        {
                            limiter::massDiffusion::markerOfTrbCellAtMatchedBC[auxUlti::getAdressOfBCEdgesOnBCValsArray(edgeId)]=false;
                        }
                    }*/
                }
                else
                    limitVal::troubleCellsMarker[nelem]=false;
            }

            totalTrbCell=parallelFuncs_Gen::sumIntOverProcs(counter);
            trouble=parallelFuncs_Gen::bool_OR_OverProcs(trouble);

            if (trouble && systemVar::currentProc==0)
            {
                std::cout << "Mass Diffusion limiter: Initial no of cells = " << totalTrbCell <<", ";
            }

            return trouble;
        }

        /**
         * @brief Function checks running condition of mass diffusion limiter using PerssonPeraireDetector.
         *
         * Basically be the same as checkTroubleCells function.
         *
         * @return
         */
        std::tuple<bool, int> checkRunning()
        {
            bool run(false);
            int numTroubleCells(0), totalTrbCell(0);
            for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
            {
                //Condition cua trouble cell la min_rho va min_rhoE phai > 0
                double minRho(limiter::massDiffusion::mathFuncs::calcMinUWOLimiter(nelem,1)), meanRhoe(0.0);
                std::tie(std::ignore,meanRhoe)=limiter::positivityPreserving::mathFuncs::calcMinMeanRhoe(nelem,1.0);

                if (minRho<0 || meanRhoe<0)
                {
                    limitVal::troubleCellsMarker[nelem]=true;
                    numTroubleCells++;
                    run=true;

                    //Check xem co trouble cell co edge la bien match hay khong, neu co thi moi active cac ham synch data
                    /*
                    int elemType(auxUlti::checkType(nelem)), edgeId(-1);
                    for (int i=0; i<elemType; i++)
                    {
                        edgeId = meshVar::inedel[nelem][i];
                        if (auxUlti::getBCType(i)==4) //Neu type = 4 thi la bien matched
                        {
                            //Mark edge cua trouble cell
                            limiter::massDiffusion::markerOfTrbCellAtMatchedBC[auxUlti::getAdressOfBCEdgesOnBCValsArray(edgeId)]=true;
                        }
                        else
                        {
                            limiter::massDiffusion::markerOfTrbCellAtMatchedBC[auxUlti::getAdressOfBCEdgesOnBCValsArray(edgeId)]=false;
                        }
                    }*/
                }
                else
                {
                    limitVal::troubleCellsMarker[nelem]=false;
                }
            }

            totalTrbCell=parallelFuncs_Gen::sumIntOverProcs(numTroubleCells);
            run=parallelFuncs_Gen::bool_OR_OverProcs(run);

            return std::make_tuple(run, totalTrbCell);
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

        /**
         * @brief Function limits Rho at only trouble cells by applying Positivity Preserving Limiter.
         *
         * This is a modified version of limitRho_PositivityPreserving in order to run at trouble cell only.
         */
        void limitRho_PPLimitiedVer()
        {
            if (mathVar::orderElem!=0)
            {
                //if (limitVal::PositivityPreserving)
                //{
                // Always run, independent of "full" version of Positivity Preserving
                    for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                    {
                        if (limitVal::troubleCellsMarker[nelem]) //condition
                            theta1Arr[nelem] = limiter::positivityPreserving::calcTheta1(nelem);
                    }
                //}
            }
        }


        namespace mathFuncs {
            void limiter_1step(int RKOrder)
            {
                //STEP 1: Calculate volume and surface rho to solve Sm = div(rho)
                limiter::massDiffusion::mathFuncs::calcVolumeGaussRho();
                limiter::massDiffusion::mathFuncs::calcRhoAtElementSurf();
                limiter::massDiffusion::mathFuncs::solveDivRho();

                //STEP 2: Calculates convective and diffusive mass flux
                limiter::massDiffusion::mathFuncs::calcSurfaceConvDiffFlux();

                //STEP 3: Solve mass conversion equation
                limiter::massDiffusion::mathFuncs::solveMassEquation(RKOrder);
            }

            void calcVolumeGaussRho()
            {
                double a(0.0), b(0.0);
                if (mathVar::orderElem>0)
                {
                    for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
                    {
                        if (limitVal::troubleCellsMarker[nelem])
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
                }
                /* Neu order = 0 thi khong can goi ham pointUVars*/
                else
                {
                    for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
                    {
                        if (limitVal::troubleCellsMarker[nelem])
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
            }

            void calcRhoAtElementSurf()
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
                            if (limitVal::troubleCellsMarker[masterCell] || limitVal::troubleCellsMarker[slaveCell])
                            {
                                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                                {
                                    std::tie(surfaceFields::rho[iedge][nG], surfaceFields::rho[iedge][nG + mathVar::nGauss1D + 1])
                                            = math::internalSurfaceValue(iedge, masterCell, nG, 1, 2);
                                }
                            }
                        }
                        else
                        {
                            std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                            if (limitVal::troubleCellsMarker[masterCell])
                            {
                                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                                {
                                    surfaceFields::rho[iedge][nG] = math::plusSideSurfaceValue(iedge, masterCell, nG, 1, 2);
                                }
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
                            if (limitVal::troubleCellsMarker[masterCell] || limitVal::troubleCellsMarker[slaveCell])
                            {
                                std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                                {
                                    //rho
                                    surfaceFields::rho[iedge][nG] = rho[masterCell][0];
                                    surfaceFields::rho[iedge][nG + mathVar::nGauss1D + 1] = rho[slaveCell][0];
                                }
                            }
                        }
                        else
                        {
                            std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                            if (limitVal::troubleCellsMarker[masterCell])
                            {
                                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                                {
                                    //rho
                                    surfaceFields::rho[iedge][nG] = rho[masterCell][0];
                                }
                            }
                        }
                    }
                }

                //Synch surface flux
                if (systemVar::parallelMode)
                    parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::rho,surfaceFields::rho);
                    //limiter_parallelFuncs::massDiff::sendRecvSurfGaussArray(surfaceFields::rho,surfaceFields::rho);
            }

            /**
             * @brief Function solves \f$S_m = \nabla (\rho)\f$ equation at trouble cells.
             */
            void solveDivRho()
            {
                std::vector<double> rhoRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                    rhoRHSTermOyDir(mathVar::orderElem + 1, 0.0);

                for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
                {
                    if (limitVal::troubleCellsMarker[nelem])
                    {
                        //2) Calculate Right hand side terms
                       limiter::massDiffusion::mathFuncs::CalcRHSTermOfDivRhoEqn(nelem, rhoRHSTermOxDir, rhoRHSTermOyDir);

                        //3) Solve for div(rho)
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            //O buoc nay, div(rho) (luu y khong co nhan mu) vao array BR1Vars::rhoX va BR1Vars::rhoY de co the dung ham pointAuxValue tinh gia tri dao ham (for convenience)
                            //Ox direction
                            BR1Vars::rhoX[nelem][iorder] = rhoRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelem][iorder];

                            //Oy direction
                            BR1Vars::rhoY[nelem][iorder] = rhoRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelem][iorder];
                        }
                    }
                }
            }

            /**
             * @brief Function calculates right hand side terms of equation \f$S_m = \nabla (\rho)\f$.
             *
             * Function is used when mass diffusion = ON.
             *
             * @param element: element Id
             * @param rhoRHSOx: right hand side term of equation \f$S_mX = \frac{\partial \rho}{\partial x}\f$. Array is returned using 'call by reference'.
             * @param rhoRHSOy: right hand side term of equation \f$S_mY = \frac{\partial \rho}{\partial y}\f$. Array is returned using 'call by reference'.
             */
            void CalcRHSTermOfDivRhoEqn(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy)
            {
                std::vector<double>
                    //vectors for volume integrals
                    rhoVolIntOx(mathVar::orderElem + 1, 0.0),
                    rhoVolIntOy(mathVar::orderElem + 1, 0.0),

                    //vectors for surface integrals
                    rhoSurfIntOx(mathVar::orderElem + 1, 0.0),
                    rhoSurfIntOy(mathVar::orderElem + 1, 0.0);

                /*1. Calculate volume integral term*/
                if (mathVar::orderElem>0)
                    limiter::massDiffusion::mathFuncs::calcVolumeIntegralTermsOfDivRhoEqn(element, rhoVolIntOx, rhoVolIntOy);

                /*2. Calculate surface integral term*/
                limiter::massDiffusion::mathFuncs::calcSurfaceIntegralTermsOfDivRhoEqn(element, rhoSurfIntOx, rhoSurfIntOy);

                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    rhoRHSOx[order] = -rhoVolIntOx[order] + rhoSurfIntOx[order];
                    rhoRHSOy[order] = -rhoVolIntOy[order] + rhoSurfIntOy[order];
                }
            }

            /**
             * @brief Function calculates volume integral term at right hand side of equation \f$S_m = \nabla (\rho)\f$.
             *
             * The term is
             * \f[-\int_{I} \nabla \phi \rho dV\f]
             *
             * @param element: element Id
             * @param rhoVolIntX: Ox volume integral term \f$-\int_{I} \frac{\partial \phi}{\partial x} \rho dV\f$
             * @param rhoVolIntY: Oy volume integral term \f$-\int_{I} \frac{\partial \phi}{\partial y} \rho dV\f$
             */
            void calcVolumeIntegralTermsOfDivRhoEqn(int element, std::vector<double> &rhoVolIntX, std::vector<double> &rhoVolIntY)
            {
                std::vector<std::vector<double>> rhoGsVol(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));

                //Calculates Gauss matrix
                //rho -------------------------------------------------------------------------------------------
                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                        rhoGsVol[na][nb]=volumeFields::rhoVolGauss[element][nanb];
                    }
                }

                for (int order = 1; order <= mathVar::orderElem; order++)
                {
                    rhoVolIntX[order] = process::volumeInte(element, rhoGsVol, order, 1);
                    rhoVolIntY[order] = process::volumeInte(element, rhoGsVol, order, 2);
                }
                rhoVolIntX[0] = 0;
                rhoVolIntY[0] = 0;
            }

            /**
             * @brief Function calculates surface integral term at right hand side of equation \f$S_m = \nabla (\rho)\f$.
             *
             * * The term is
             * \f[\int_{\partial I} \phi \rho \cdot \mathbf{n} dV\f]
             *
             * @param element: element Id
             * @param rhoSurfIntX: Ox surface integral term \f$\int_{\partial I} \phi \rho n_x dV\f$
             * @param rhoSurfIntY: OY surface integral term \f$\int_{\partial I} \phi \rho n_y dV\f$
             */
            void calcSurfaceIntegralTermsOfDivRhoEqn(int element, std::vector<double> &rhoSurfIntX, std::vector<double> &rhoSurfIntY)
            {
                int elemType(auxUlti::checkType(element)), edgeName(0);
                std::vector<std::vector<double>> rhoFluxX(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0)),
                    rhoFluxY(mathVar::nGauss1D + 1, std::vector<double>(elemType, 0.0));

                std::vector<double> rhoFluxXTemp(mathVar::nGauss1D + 1, 0.0),
                    rhoFluxYTemp(mathVar::nGauss1D + 1, 0.0);

                /*1. Calculate flux of rho at all Gauss points on all faces of element*/
                limiter::massDiffusion::mathFuncs::getGaussVectorOfRho(element, rhoFluxX, rhoFluxY);

                /*2. Calculates surface integrals of rho at all order*/
                for (int nface = 0; nface < elemType; nface++)
                {
                    edgeName = meshVar::inedel[element][nface];
                    for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                    {
                        rhoFluxXTemp[nG] = rhoFluxX[nG][nface];
                        rhoFluxYTemp[nG] = rhoFluxY[nG][nface];
                    }

                    for (int order = 0; order <= mathVar::orderElem; order++)
                    {
                        rhoSurfIntX[order] += process::surfaceInte(element, edgeName, rhoFluxXTemp, order);
                        rhoSurfIntY[order] += process::surfaceInte(element, edgeName, rhoFluxYTemp, order);
                    }
                }
            }

            /**
             * @brief Function gets values  all Rho fluxes at all Gauss points on all edges of cell.
             *
             * Rho flux is calculated by central flux:
             * \f[Flux_\rho = \frac{1}{2}(\rho_+ + \rho_-) \cdot \mathbf{n}\f]
             * Function returns Rho fluxes through 'call by reference'.
             *
             * @param element: element Id
             * @param rhoFluxX: term \f${Flux_\rho}x = \frac{1}{2}(\rho_+ + \rho_-) nx\f$
             * @param rhoFluxY: term \f${Flux_\rho}y = \frac{1}{2}(\rho_+ + \rho_-) ny\f$
             */
            void getGaussVectorOfRho(int element, std::vector<std::vector<double>> &rhoFluxX, std::vector<std::vector<double>> &rhoFluxY)
            {
                int elemType(auxUlti::checkType(element)), edgeId(0);
                int faceBcType(0);
                std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0));

                for (int nface = 0; nface < elemType; nface++)
                {
                    edgeId = meshVar::inedel[element][nface];
                    faceBcType = auxUlti::getBCType(edgeId);
                    double nx(auxUlti::getNormVectorComp(element, edgeId, 1)), ny(auxUlti::getNormVectorComp(element, edgeId, 2)), rhoP(0.0), rhoM(0.0);
                    if (faceBcType == 0)  //internal edge
                    {
                        for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
                        {
                            std::tie(rhoP, rhoM) = auxUlti::getUAtInterfaces(edgeId, element, nGauss, 1);
                            rhoFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(rhoM, rhoP, nx);
                            rhoFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(rhoM, rhoP, ny);
                        }
                    }
                    else  //boundary edge
                    {
                        for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
                        {
                            std::tie(rhoP, rhoM)=rhoBCsImplement(edgeId, nGauss);
                            rhoFluxX[nGauss][nface]=math::numericalFluxes::auxFlux(rhoP, rhoM, nx);
                            rhoFluxY[nGauss][nface]=math::numericalFluxes::auxFlux(rhoP, rhoM, ny);
                        }
                    }
                }
            }

            void calcSurfaceConvDiffFlux()
            {
                //NOTE: trong ham limiter, gia tri luu trong array BR1Vars::rhoX va BR1Vars::rhoY (duoc dung trong ham pointAuxValue) la gia tri div(rho) khong phai mu*div(rho)

                //Tinh flux qua cac internal edge.
                int masterCell(-1), slaveCell(-1), bcGrp(0);
                double TMaster(0.0), TSlave(0.0), muMaster(0.0), muSlave(0.0);

                std::vector<double> UMaster(3, 0.0), USlave(3, 0.0); //only 3 rows because of not including rhoE

                double drhoXMaster(0.0), drhoYMaster(0.0), drhoXSlave(0.0), drhoYSlave(0.0), nx, ny;

                double
                        //Local fluxes
                        convectiveFluxXMaster(0.0),
                        convectiveFluxYMaster(0.0),
                        diffusiveFluxXMaster(0.0),
                        diffusiveFluxYMaster(0.0),

                        convectiveFluxXSlave(0.0),
                        convectiveFluxYSlave(0.0),
                        diffusiveFluxXSlave(0.0),
                        diffusiveFluxYSlave(0.0);

                for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
                {
                    bcGrp = auxUlti::getBCType(iedge);
                    if (bcGrp == 0)
                    {
                        std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                        if (limitVal::troubleCellsMarker[masterCell] || limitVal::troubleCellsMarker[slaveCell])
                        {
                            nx = auxUlti::getNormVectorComp(masterCell, iedge, 1);
                            ny = auxUlti::getNormVectorComp(masterCell, iedge, 2);

                            //Luon dung LxF flux cho mass diffusion limiter -> Tinh he so C cua LxF flux
                            math::numericalFluxes::findMaxLxFConstantOnEdge(iedge,masterCell);

                            for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                            {
                                std::tie(TMaster,TSlave)=auxUlti::getTAtInterfaces(iedge,masterCell,nG);

                                for (int i = 0; i < 3; i++) //Khong lay gia tri rhoE
                                {
                                    //std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
                                    std::tie(UMaster[i], USlave[i]) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, i + 1);
                                }
                                if (flowProperties::viscous)
                                {
                                    std::tie(drhoXMaster, drhoXSlave) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, 1, 1);
                                    std::tie(drhoYMaster, drhoYSlave) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, 1, 2);

                                    //Save derivaties to surfaceFields
                                    surfaceFields::dRhoX[iedge][nG]=drhoXMaster;
                                    surfaceFields::dRhoX[iedge][mathVar::nGauss1D+nG+1]=drhoXSlave;
                                    surfaceFields::dRhoY[iedge][nG]=drhoYMaster;
                                    surfaceFields::dRhoY[iedge][mathVar::nGauss1D+nG+1]=drhoYSlave;
                                }

                                //Trong limiter, dung LxF flux de do stability cao nhat
                                //Do do, khong can chuyen he toa do tu global ve local
                                convectiveFluxXMaster=UMaster[1]; //rhou
                                convectiveFluxXSlave=USlave[1]; //rhou

                                convectiveFluxYMaster=UMaster[2]; //rhov
                                convectiveFluxYSlave=USlave[2]; //rhov

                                if (flowProperties::viscous)
                                {
                                    //Phai tinh mu va lay mu*div(rho) vi khi dung mass diffusion limiter, khong giai bien phu S=mu*div(U) ma chi giai Sm=div(rho)
                                    muMaster=math::CalcVisCoef(TMaster);
                                    muSlave=math::CalcVisCoef(TSlave);
                                    diffusiveFluxXMaster=-limiter::massDiffusion::DmCoeff*muMaster*drhoXMaster/UMaster[0];
                                    diffusiveFluxYMaster=-limiter::massDiffusion::DmCoeff*muMaster*drhoYMaster/UMaster[0];
                                    diffusiveFluxXSlave=-limiter::massDiffusion::DmCoeff*muSlave*drhoXSlave/USlave[0];
                                    diffusiveFluxYSlave=-limiter::massDiffusion::DmCoeff*muSlave*drhoYSlave/USlave[0];
                                }

                                //Tinh flux
                                //Vi dang tinh cho cell master cua edge nen phia trong (phia +) la cua cell master
                                surfaceFields::invis_rho[iedge][nG]=0.5*(nx*(convectiveFluxXMaster+convectiveFluxXSlave)
                                                                         + ny*(convectiveFluxYMaster+convectiveFluxYSlave)
                                                                         - //chu y dau -
                                                                         LxFConst[iedge]*(USlave[0]-UMaster[0]));
                                surfaceFields::invis_rho[iedge][nG + mathVar::nGauss1D + 1]=-surfaceFields::invis_rho[iedge][nG];

                                //Viscous flux tinh bang central flux
                                if (flowProperties::viscous)
                                {
                                    surfaceFields::Vis_rho[iedge][nG]=0.5*(nx*(diffusiveFluxXMaster+diffusiveFluxXSlave)
                                                                               + ny*(diffusiveFluxYMaster+diffusiveFluxYSlave));
                                    surfaceFields::Vis_rho[iedge][nG + mathVar::nGauss1D + 1]=-surfaceFields::Vis_rho[iedge][nG];
                                }
                            }
                        }
                    }
                    else
                    {
                        if (flowProperties::viscous)
                        {
                            //int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(iedge));
                            std::tie(masterCell, std::ignore) = auxUlti::getMasterServantOfEdge(iedge);
                            if (limitVal::troubleCellsMarker[masterCell] || limitVal::troubleCellsMarker[slaveCell])
                            {
                                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                                {
                                    std::tie(TMaster,std::ignore)=auxUlti::getTAtInterfaces(iedge,masterCell,nG);
                                    muMaster=math::CalcVisCoef(TMaster);
                                    surfaceFields::dRhoX[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,1,1)*muMaster;
                                    surfaceFields::dRhoY[iedge][nG] = math::plusSideSurfaceDerivativeValue(iedge,masterCell,nG,1,2)*muMaster;
                                }
                            }
                        }
                    }
                }

                //Synch surface flux
                if (systemVar::parallelMode && flowProperties::viscous)
                {
                    parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhoX, surfaceFields::dRhoX);
                    parallelFuncs_GaussPt::sendRecvSurfGaussArray(surfaceFields::dRhoY, surfaceFields::dRhoY);
                    //limiter_parallelFuncs::massDiff::sendRecvSurfGaussArray(surfaceFields::dRhoX, surfaceFields::dRhoX);
                    //limiter_parallelFuncs::massDiff::sendRecvSurfGaussArray(surfaceFields::dRhoY, surfaceFields::dRhoY);
                }
            }

            std::tuple<double, double, double, double> calcGaussConvDiffTerms(int element, int na, int nb)
            {
                double convTermX,diffTermX,convTermY,diffTermY;

                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                std::vector<std::vector<double>> InviscidTerm(4, std::vector<double>(2, 0.0));
                double
                    rhoVal(volumeFields::rhoVolGauss[element][nanb]),
                    rhouVal(volumeFields::rhouVolGauss[element][nanb]),
                    rhovVal(volumeFields::rhovVolGauss[element][nanb]),
                    a(0.0),b(0.0),muVal(0.0),dRhoX,dRhoY;

                std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                muVal=math::CalcVisCoef(volumeFields::T[element][nanb]);

                dRhoX=math::pointAuxValue(element,a,b,1,1);
                dRhoY=math::pointAuxValue(element,a,b,1,2);

                convTermX=rhouVal;
                convTermY=rhovVal;
                diffTermX=-limiter::massDiffusion::DmCoeff*muVal*dRhoX/rhoVal; //Phai nhan mu vi bien phu luc nay la div(rho)
                diffTermY=-limiter::massDiffusion::DmCoeff*muVal*dRhoY/rhoVal; //Phai nhan mu vi bien phu luc nay la div(rho)

                return std::make_tuple(convTermX, diffTermX, convTermY, diffTermY);
            }

            void calcVolumeIntegralConvDiffTerms(int element, std::vector<double>&VolIntTerm1)
            {
                std::vector<std::vector<double>>
                    GsVolX1(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
                    GsVolY1(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));
                double invTermX(0.0), visTermX(0.0), invTermY(0.0), visTermY(0.0);

                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        std::tie(invTermX,visTermX,invTermY,visTermY)=limiter::massDiffusion::mathFuncs::calcGaussConvDiffTerms(element,na,nb);

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
                std::vector<std::vector<double>>massFlux(mathVar::nGauss1D + 1, std::vector<double>(4, 0.0));
                std::vector<double> massFluxBC(4, 0.0), massFluxTemp(mathVar::nGauss1D + 1, 0.0);

                double convTerm(0.0), diffTerm(0.0);

                //std::vector<double> SurInt(4, 0.0);

                for (int nface = 0; nface < elemType; nface++)
                {
                    edgeName = meshVar::inedel[element][nface];
                    faceBcType = auxUlti::getBCType(edgeName);

                    if (faceBcType == 0)  //internal edge
                    {
                        for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
                        {
                            /*INVISCID TERM*/
                            //Get value
                            std::tie(convTerm, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 1, 1);

                            /*VISCOUS TERM*/
                            //Get value
                            if (flowProperties::viscous)
                            {
                                std::tie(diffTerm, std::ignore) = process::NSFEq::getInterfacesFluxes(edgeName, element, nGauss, 2, 1);
                            }

                            /*Calculate fluxes*/
                            massFlux[nGauss][nface] = convTerm+diffTerm;
                        }
                    }
                    else  //boundary edge
                    {
                        for (int nGauss = 0; nGauss <= mathVar::nGauss1D; nGauss++)
                        {
                            massFluxBC = NSFEqBCsImplement(element, edgeName, nGauss);
                            massFlux[nGauss][nface] = massFluxBC[0];
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
                            massFluxTemp[nG] = massFlux[nG][nface];
                        }
                        SurfIntTerm1[order] += process::surfaceInte(element, edgeName, massFluxTemp, order);
                    }
                }
            }

            void solveMassEquation(int RKOrder)
            {
                std::vector<double>
                    RHSTerm1(mathVar::orderElem + 1, 0.0),
                    ddtRhoVector(mathVar::orderElem + 1, 0.0),
                    rhoVectorN(mathVar::orderElem + 1, 0.0),
                    UnVector(mathVar::orderElem + 1, 0.0);

                for (int nelem=0; nelem<meshVar::nelem2D; nelem++)
                {
                    if (limitVal::troubleCellsMarker[nelem])
                    {
                        //2) Calculate Right hand side terms
                        std::vector<double>
                            VolIntTerm1(mathVar::orderElem + 1, 0.0),
                            SurfIntTerm1(mathVar::orderElem + 1, 0.0);

                        /*Volume integral term===========================================================================*/
                        limiter::massDiffusion::mathFuncs::calcVolumeIntegralConvDiffTerms(nelem, VolIntTerm1);

                        /*Surface integral term===========================================================================*/
                        limiter::massDiffusion::mathFuncs::calcSurfaceIntegralTerms(nelem, SurfIntTerm1);

                        for (int order = 0; order <= mathVar::orderElem; order++)
                        {
                            RHSTerm1[order] = VolIntTerm1[order] - SurfIntTerm1[order];
                        }

                        //3) Solve for time derivartives of conservative variables
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            ddtRhoVector[iorder] = RHSTerm1[iorder] / stiffMatrixCoeffs[nelem][iorder];
                        }

                        //4) Solve time marching
                        //rho
                        for (int order = 0; order <= mathVar::orderElem; order++)
                        {
                            UnVector[order] = rho[nelem][order];
                        }
                        rhoVectorN = process::NSFEq::solveTimeMarching(nelem, ddtRhoVector, UnVector, RKOrder, 1);

                        //5) Save results to conservative variables arrays
                        for (int order = 0; order <= mathVar::orderElem; order++)
                        {
                            rhoN[nelem][order] = rhoVectorN[order];
                        }
                    }
                }
            }

            void calcVolIntegralTerm(int element,std::vector<double>&VolIntTerm1)
            {
                std::vector<std::vector<double>> ViscousTerms(4, std::vector<double>(2, 0.0));
                std::vector<std::vector<double>> InviscidTerms(4, std::vector<double>(2, 0.0));
                //std::vector<double> VolInt(4, 0.0);

                std::vector<std::vector<double>>
                    GsVolX1(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
                    GsVolY1(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));

                double rhoVal,rhouVal,rhovVal,uVal(0.0),vVal(0.0),a(0.0),b(0.0),muVal(0.0),dRhoX(0.0),dRhoY(0.0);

                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));

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
                        GsVolX1[na][nb] += -limiter::massDiffusion::DmCoeff*dRhoX/rhoVal;
                        /*B2. Viscous term on Oy direction*/
                        GsVolY1[na][nb] += -limiter::massDiffusion::DmCoeff*dRhoY/rhoVal;
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

            double calcMinRhoResidual(int element)
            {
                //std::vector<double> vectorRho((mathVar::nGauss+1)*(mathVar::nGauss+5), 1.0),
                        //vectorRho0((mathVar::nGauss+1)*(mathVar::nGauss+5), 1.0);
                std::vector<double> vectorRho((mathVar::nGauss2D+1)*(mathVar::nGauss2D+1)+(mathVar::nGauss1D+1)*4, 1.0),
                        vectorRho0((mathVar::nGauss2D+1)*(mathVar::nGauss2D+1)+(mathVar::nGauss1D+1)*4, 1.0);

                double aG(0.0), bG(0.0), minRho(0.0), minRho0(0.0);
                int counter(0);

                //Compute rho at all internal Gauss point

                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);
                        vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                        vectorRho0[counter]=(limiter::massDiffusion::mathFuncs::pointRho0(element, aG, bG));
                        counter++;
                    }
                }

                //Compute rho at edge DA
                aG = -1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(limiter::massDiffusion::mathFuncs::pointRho0(element, aG, bG));
                    counter++;
                }
                //Compute rho at edge BC
                aG = 1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(limiter::massDiffusion::mathFuncs::pointRho0(element, aG, bG));
                    counter++;
                }
                //Compute rho at edge AB
                bG = -1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(limiter::massDiffusion::mathFuncs::pointRho0(element, aG, bG));
                    counter++;
                }
                //Compute rho at edge CD
                bG = 1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    vectorRho[counter]=(math::pointValueNoLimiter(element, aG, bG, 1));
                    vectorRho0[counter]=(limiter::massDiffusion::mathFuncs::pointRho0(element, aG, bG));
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

            double calcMinUWOLimiter(int element, int type)
            {
                double aG(0.0), bG(0.0), min(1e10), val(0.0);
                int counter(0);

                //Compute rho at all internal Gauss point

                for (int na = 0; na <= mathVar::nGauss2D; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                    {
                        std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);
                        val=(math::pointValueNoLimiter(element, aG, bG, type));
                        counter++;
                        if (val<min)
                            min=val;
                    }
                }

                //Compute rho at edge DA
                aG = -1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    val=(math::pointValueNoLimiter(element, aG, bG, type));
                    counter++;
                    if (val<min)
                        min=val;
                }
                //Compute rho at edge BC
                aG = 1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    val=(math::pointValueNoLimiter(element, aG, bG, type));
                    counter++;
                    if (val<min)
                        min=val;
                }
                //Compute rho at edge AB
                bG = -1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    val=(math::pointValueNoLimiter(element, aG, bG, type));
                    counter++;
                    if (val<min)
                        min=val;
                }
                //Compute rho at edge CD
                bG = 1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    val=(math::pointValueNoLimiter(element, aG, bG, type));
                    if (val<min)
                        min=val;
                    counter++;
                }
                return min;
            }

            double calcMaxAbsGradRhoOnEdges(int element)
            {
                double aG(0.0), bG(0.0), max(1e-15), val(0.0);

                //Compute rho at edge DA
                aG = -1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    bG = mathVar::xGaussSur[nG];

                    //Tinh do lon cua grad(rho)
                    val=pow(pow(math::pointAuxValue(element, aG, bG, 1, 1),2) + pow(math::pointAuxValue(element, aG, bG, 1, 2),2),0.5);
                    if (val>max)
                        max=val;
                }
                //Compute rho at edge BC
                aG = 1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    val=pow(pow(math::pointAuxValue(element, aG, bG, 1, 1),2) + pow(math::pointAuxValue(element, aG, bG, 1, 2),2),0.5);
                    if (val>max)
                        max=val;
                }
                //Compute rho at edge AB
                bG = -1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    val=pow(pow(math::pointAuxValue(element, aG, bG, 1, 1),2) + pow(math::pointAuxValue(element, aG, bG, 1, 2),2),0.5);
                    if (val>max)
                        max=val;
                }
                //Compute rho at edge CD
                bG = 1;
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    val=pow(pow(math::pointAuxValue(element, aG, bG, 1, 1),2) + pow(math::pointAuxValue(element, aG, bG, 1, 2),2),0.5);
                    if (val>max)
                        max=val;
                }
                return max;
            }

            double calcMeanAbsGradRho(int element)
            {
                return
                        (pow(
                             pow(BR1Vars::rhoX[element][0],2)+
                             pow(BR1Vars::rhoY[element][0],2)
                             ,0.5));
            }
        }
    }

    namespace IOMassDiff {
        void readSetting()
        {
            //Input chieu dai cac array chua data can doc o day, luu y loai dataType nao khong can doc, van phai de chieu day array la 1
            const int doubleArrSize(2);
            const int intArrSize(1);
            const int boolArrSize(1);
            const int strArrSize(1);

            //Intput so luong cac bien can doc o day
            int numDbl(2),
                    numInt(1),
                    numBool(0),
                    numStr(0);

            std::string file("LimiterSettings.txt");
            std::string loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
            std::string
                    keyWordsDouble[doubleArrSize] = { "Co", "DmCoeff" },
                    keyWordsInt[intArrSize] = {"maxIter"},
                    keyWordsBool[boolArrSize] = {},
                    keyWordsStr[strArrSize] = {};

            double outDB[doubleArrSize] = {};
            int outInt[intArrSize] = {};
            bool outBool[boolArrSize] = {};
            std::string outStr[strArrSize] = {};

            IO::readDataFile(file, loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr,
                             outDB, outInt, outBool, outStr,
                             numDbl, numInt, numBool, numStr);

            //Save data to variables
            limiter::massDiffusion::pseudoCo=outDB[0];
            limiter::massDiffusion::DmCoeff=outDB[1];
            limiter::massDiffusion::maxIter=outInt[0];
        }
    }
}
