#include "positivityPreserving.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include <mpi.h>
#include <math.h>
#include <tuple>
#include <iostream>
#include <vector>
#include <algorithm>

namespace limiter
{
    namespace positivityPreserving
    {
        void limiter()
        {
            if (limitVal::PositivityPreserving)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {
                    if (!flowProperties::massDiffusion)
                    {
                        theta1Arr[nelem] = limiter::positivityPreserving::calcTheta1(nelem);
                    }
                    theta2Arr[nelem] = limiter::positivityPreserving::calcTheta2(nelem);
                }

                if (systemVar::currentProc==0)
                {
                    int numOfTotalLimitedCell(limitVal::numOfLimitCell), recvNum;
                    for (int irank=1;irank<systemVar::totalProc;irank++) {
                        MPI_Recv(&recvNum, 1, MPI_INT, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        numOfTotalLimitedCell+=recvNum;
                    }
                    if (numOfTotalLimitedCell > 0)
                    {
                        std::cout << "Posivity preserving limiter is applied at " << numOfTotalLimitedCell << " cell(s)\n";
                    }
                }
                else {
                    MPI_Send(&limitVal::numOfLimitCell, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
                limitVal::numOfLimitCell = 0;
            }
        }

        void initialiseThetaVector()
        {
            for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
            {
                theta1Arr[nelem] = 1.0;
                theta2Arr[nelem] = 1.0;
            }
        }

        double calcTheta2(int element)
        {
            //Find theta2
            double meanRho(0.0), theta1(theta1Arr[element]), theta2(0.0), meanRhoe(0.0), minRhoe(0.0);

            std::tie(minRhoe, meanRhoe) = limiter::positivityPreserving::mathFuncs::calcMinMeanRhoe(element, theta1);
            theta2 = limiter::positivityPreserving::mathFuncs::calcTheta2Coeff(meanRhoe, minRhoe, meanRho);

            if ((theta2 < 1))
            {
                limitVal::numOfLimitCell++;
            }
            //Reset limit flag
            limitVal::limitFlagLocal = false;
            return theta2;
        }

        double calcTheta1(int element)
        {
            double meanRho(0.0), minRho(0.0), theta1(0.0);

            //Find theta1
            minRho = limiter::positivityPreserving::mathFuncs::calcMinRho(element);
            meanRho = rho[element][0];

            //Compute theta1
            theta1 = limiter::positivityPreserving::mathFuncs::calcTheta1Coeff(minRho, meanRho);
            return theta1;
        }

        namespace mathFuncs
        {
            double calcMinRho(int element)
            {
                double aG(0.0), bG(0.0), min(1e10), rhoVal(0.0);
                int counter(0);

                //Compute rho at all internal Gauss point

                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);
                        rhoVal=(math::pointValueNoLimiter(element, aG, bG, 1));
                        counter++;
                        if (rhoVal<min)
                            min=rhoVal;
                    }
                }

                //Compute rho at edge DA
                aG = -1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    rhoVal=(math::pointValueNoLimiter(element, aG, bG, 1));
                    counter++;
                    if (rhoVal<min)
                        min=rhoVal;
                }
                //Compute rho at edge BC
                aG = 1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    bG = mathVar::xGaussSur[nG];
                    rhoVal=(math::pointValueNoLimiter(element, aG, bG, 1));
                    counter++;
                    if (rhoVal<min)
                        min=rhoVal;
                }
                //Compute rho at edge AB
                bG = -1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    rhoVal=(math::pointValueNoLimiter(element, aG, bG, 1));
                    counter++;
                    if (rhoVal<min)
                        min=rhoVal;
                }
                //Compute rho at edge CD
                bG = 1;
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    aG = mathVar::xGaussSur[nG];
                    rhoVal=(math::pointValueNoLimiter(element, aG, bG, 1));
                    if (rhoVal<min)
                        min=rhoVal;
                    counter++;
                }
                return min;
            }

            //Function calculates modified value of Rho at abitrary point (for calculating theta2)
            double calcRhoModified(int element, double a, double b, double theta1)
            {
                double rhoMod(0.0);
                std::vector<double> Value(mathVar::orderElem + 1, 0.0);

                Value = auxUlti::getElementConserValuesOfOrder(element, 1);
                math::basisFc(a, b, auxUlti::checkType(element));
                for (int order = 1; order <= mathVar::orderElem; order++)
                {
                    rhoMod += Value[order] * mathVar::B[order] * theta1;
                }
                rhoMod += Value[0];
                return rhoMod;
            }

            std::tuple<double, double> calcMinMeanRhoe(int element, double theta1)
            {
                /* Neu co mass diffusion, phai modify lai rhou va rhov thanh rhou_m va rhov_m.
                 * Vi khi giai T explicit, ham CalcTFromConsvVar_massDiff tinh gia tri mass
                 * diffusion flux Jd = (1/Sc)*mu*div(rho)/rho voi Sc la so Schmidt (gia tri DmCoeff
                 * la 1/Sc) va mu tinh tu T cua step truoc (T_old), nen de limit dung can phai tinh
                 * Jd o ham limiter giong nhu o ham CalcTFromConsvVar_massDiff.
                */

                double aG(0.0), bG(0.0), min(1e10), rhoMod(0.0), rhouOrg(0.0), rhovOrg(0.0), rhoEOrg(0.0), sum(0), mean(0), rhoeVal(0.0);
                int counter(0), edgeId;

                //Compute rho at all internal Gauss point
                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            int nanb(calcArrId(na,nb,mathVar::nGauss+1));
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(volumeFields::T[element][nanb]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod,rhouOrg,rhovOrg,rhoEOrg);
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                }

                if (auxUlti::checkType(element)==4) // cell tu giac
                {
                    //Tinh rhoe tren edge
                    int edgeAB(0), edgeBC(1), edgeCD(2), edgeDA(3);

                    //Compute rho at edge DA
                    aG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeDA);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=(limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                    //Compute rho at edge BC
                    aG = 1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeBC);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=(limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                    //Compute rho at edge AB
                    bG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        aG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeAB);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=(limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                    //Compute rho at edge CD

                    bG = 1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        aG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeCD);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=(limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                }
                else if (auxUlti::checkType(element)==3) // cell tam giac
                {
                    //Tinh rhoe tren edge
                    int edgeAB(0), edgeBC(1), edgeCA(2);

                    //Compute rho at edge CA
                    aG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeCA);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=(limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                    //Compute rho at edge BC
                    aG = 1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        bG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeBC);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=(limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                    //Compute rho at edge AB
                    bG = -1;
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        aG = mathVar::xGaussSur[nG];
                        rhoMod = limiter::positivityPreserving::mathFuncs::calcRhoModified(element, aG, bG, theta1);
                        rhouOrg = math::pointValueNoLimiter(element, aG, bG, 2);
                        rhovOrg = math::pointValueNoLimiter(element, aG, bG, 3);
                        rhoEOrg = math::pointValueNoLimiter(element, aG, bG, 4);

                        if (flowProperties::massDiffusion)
                        {
                            edgeId=auxUlti::getEdgeHasInputOrderOfElement(element,edgeAB);
                            double rhoX(math::pointAuxValue(element,aG,bG,5,1)), rhoY(math::pointAuxValue(element,aG,bG,5,2)),
                                    mu(math::CalcVisCoef(surfaceFields::T[edgeId][nG]));
                            rhouOrg=rhouOrg-material::massDiffusion::DmCoeff*rhoX*mu/rhoMod;
                            rhovOrg=rhovOrg-material::massDiffusion::DmCoeff*rhoY*mu/rhoMod;
                        }
                        rhoeVal=(limiter::positivityPreserving::mathFuncs::calRhoeFromConserVars(rhoMod, rhouOrg, rhovOrg, rhoEOrg));
                        if (rhoeVal<min)
                            min=rhoeVal;
                        sum+=rhoeVal;
                        counter++;
                    }
                }

                mean = sum/counter;
                return std::make_tuple(min, mean);
            }

            double calRhoeFromConserVars(double rho, double rhou, double rhov, double rhoE)
            {
                double rhoe(0.0);
                rhoe = rhoE - 0.5*(rhou*rhou + rhov * rhov) / rho;
                return rhoe;
            }

            double calcTheta1Coeff(double minRho, double meanRho)
            {
                double temp1(0.0), theta(0.0);
                std::vector<double> vectorOmega(2, 0.0);
                vectorOmega[0] = systemVar::epsilon;
                vectorOmega[1] = meanRho;
                double omega(*std::min_element(vectorOmega.begin(), vectorOmega.end()));  //find min value of vector

                //Compare meanRho and minRho, if they are equal then theta = 1.0 (means no limiter applied)
                if (fabs(meanRho - minRho/meanRho) <= 1e-8)
                {
                    theta = 1.0;
                }
                else
                {
                    temp1 = (meanRho - omega) / (meanRho - minRho);
                    if (temp1 < 1.0)
                    {
                        theta = temp1;
                    }
                    else
                    {
                        theta = 1.0;
                    }
                }
                return theta;
            }

            double calcTheta2Coeff(double meanRhoe, double minRhoe, double meanRho)
            {
                double temp1(0.0), theta(0.0);
                std::vector<double> vectorOmega(3, 0.0);
                vectorOmega[0] = systemVar::epsilon;
                vectorOmega[1] = meanRho;
                vectorOmega[2] = meanRhoe;
                double omega(*std::min_element(vectorOmega.begin(), vectorOmega.end()));  //find min value of vector

                //Compare meanRhoe and minRhoe, if they are equal then theta = 1.0 (means no limiter applied)
                if (fabs(meanRhoe - minRhoe/meanRhoe) <= 1e-8)
                {
                    theta = 1.0;
                }
                else
                {
                    temp1 = (meanRhoe - omega) / (meanRhoe - minRhoe);
                    if (temp1 < 1.0)
                    {
                        theta = temp1;
                    }
                    else
                    {
                        theta = 1.0;
                    }
                }

                return theta;
            }
        }
    }
}
