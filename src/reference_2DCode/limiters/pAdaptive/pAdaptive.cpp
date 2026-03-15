#include "pAdaptive.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include "./limiters/mathFunctions.h"
#include "DGMath.h"
#include <math.h>

namespace limiter
{
    namespace pAdaptive
    {
        void limiter()
        {
            //p-Adaptive
            if (limitVal::PAdaptive)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {

                    for (int valType = 1; valType <= 4; valType++)
                    {
                        //p-Adaptive limiter
                        limiter::pAdaptive::limiter_1Elem(nelem, valType);
                    }
                    if (limitVal::pAdaptive::limitFlagLocal == true)
                    {
                        limitVal::pAdaptive::numOfLimitCell++;
                        limitVal::pAdaptive::limitFlagLocal = false;
                    }
                }
                if (limitVal::pAdaptive::numOfLimitCell > 0)
                {
                    std::cout << "P-adaptive limiter is applied at " << limitVal::pAdaptive::numOfLimitCell << " cell(s)\n";
                    limitVal::pAdaptive::numOfLimitCell = 0;
                }
            }
        }

        void limiter_1Elem(int element, int valType)
        {
            int elemType(auxUlti::checkType(element)),
                elemJPlus(0), elemJMinus(0),  //means j+1, j-1
                elemIPlus(0), elemIMinus(0);  //means i+1, i-1
            switch (elemType)
            {
            case 3:
            {
                std::vector<double> outputVar(mathVar::orderElem, 0.0);
                elemIPlus = meshVar::neighboringElements[element][1];
                elemIMinus = meshVar::neighboringElements[element][2];
                if (mathVar::orderElem > 1)
                {
                    elemJMinus = meshVar::neighboringElements[element][0];
                }

                if ((elemIPlus >= 0) & (elemIMinus >= 0) & (elemJMinus >= 0))
                {
                    switch (valType)
                    {
                    case 1:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rho[elemIPlus][0], rho[elemIMinus][0], rho[elemJMinus][0]);
                        rho[element][1] = outputVar[0];
                        rho[element][2] = outputVar[1];
                    }
                    break;
                    case 2:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rhou[elemIPlus][0], rhou[elemIMinus][0], rhou[elemJMinus][0]);
                        rhou[element][1] = outputVar[0];
                        rhou[element][2] = outputVar[1];
                    }
                    break;
                    case 3:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rhov[elemIPlus][0], rhov[elemIMinus][0], rhov[elemJMinus][0]);
                        rhov[element][1] = outputVar[0];
                        rhov[element][2] = outputVar[1];
                    }
                    break;
                    case 4:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Tri(element, valType, rhoE[elemIPlus][0], rhoE[elemIMinus][0], rhoE[elemJMinus][0]);
                        rhoE[element][1] = outputVar[0];
                        rhoE[element][2] = outputVar[1];
                    }
                    break;
                    default:
                        break;
                    }
                }
            }
            break;
            case 4:
            {
                std::vector<double> outputVar(mathVar::orderElem, 0.0);
                //bool internalElem(true);
                elemIPlus = meshVar::neighboringElements[element][1];
                elemIMinus = meshVar::neighboringElements[element][3];
                if (mathVar::orderElem > 1)
                {
                    elemJPlus = meshVar::neighboringElements[element][2];
                    elemJMinus = meshVar::neighboringElements[element][0];
                }

                if ((elemIPlus >= 0) & (elemIMinus >= 0) & (elemJPlus >= 0) & (elemJMinus >= 0))
                {
                    switch (valType)
                    {
                    case 1:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rho[elemIPlus][0], rho[elemIMinus][0], rho[elemJPlus][0], rho[elemJMinus][0]);
                        rho[element][1] = outputVar[0];
                        rho[element][2] = outputVar[1];
                    }
                    break;
                    case 2:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rhou[elemIPlus][0], rhou[elemIMinus][0], rhou[elemJPlus][0], rhou[elemJMinus][0]);
                        rhou[element][1] = outputVar[0];
                        rhou[element][2] = outputVar[1];
                    }
                    break;
                    case 3:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rhov[elemIPlus][0], rhov[elemIMinus][0], rhov[elemJPlus][0], rhov[elemJMinus][0]);
                        rhov[element][1] = outputVar[0];
                        rhov[element][2] = outputVar[1];
                    }
                    break;
                    case 4:
                    {
                        outputVar = limiter::pAdaptive::pAdaptiveChildFunction_Quad(element, valType, rhoE[elemIPlus][0], rhoE[elemIMinus][0], rhoE[elemJPlus][0], rhoE[elemJMinus][0]);
                        rhoE[element][1] = outputVar[0];
                        rhoE[element][2] = outputVar[1];
                    }
                    break;
                    default:
                        break;
                    }
                }
            }
            break;
            default:
                break;
            }
        }

        std::vector<double> pAdaptiveChildFunction_Quad(int element, int valType, double IPlus, double IMinus, double JPlus, double JMinus)
        {
            std::vector<double>elementConsValOfOrder(mathVar::orderElem + 1, 0.0),
                inputArgument(3, 0.0), output(mathVar::orderElem, 0.0);
            elementConsValOfOrder = auxUlti::getElementConserValuesOfOrder(element, valType);

            double UBC(math::pointValueNoLimiter(element, 1.0, 0.0, valType) - elementConsValOfOrder[0]),
                UAD(math::pointValueNoLimiter(element, -1.0, 0.0, valType) - elementConsValOfOrder[0]),
                UCD(math::pointValueNoLimiter(element, 0.0, 1.0, valType) - elementConsValOfOrder[0]),
                UAB(math::pointValueNoLimiter(element, 0.0, -1.0, valType) - elementConsValOfOrder[0]),
                UBCMod(0.0), UADMod(0.0), UCDMod(0.0), UABMod(0.0),
                UBC_check(0.0), UAD_check(0.0), UCD_check(0.0), UAB_check(0.0),
                M(limiter::mathForLimiter::calM(element, valType)), Lxy(meshVar::localCellSize[element]);

            inputArgument[0] = UBC;
            inputArgument[1] = IPlus - elementConsValOfOrder[0];
            inputArgument[2] = elementConsValOfOrder[0] - IMinus;
            UBCMod = limiter::mathForLimiter::minmod(inputArgument);
            UBC_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1); //*fabs(elementConsValOfOrder[0])

            inputArgument[0] = -UAD;
            UADMod = -limiter::mathForLimiter::minmod(inputArgument);
            inputArgument[0] = UAD;
            UAD_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);  //*fabs(elementConsValOfOrder[0])
            if (fabs(UBC - UBC_check)>1e-10 || fabs(UAD - UAD_check)>1e-10)
            {
                //std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
                if (limitVal::pAdaptive::limitFlagLocal == false)
                {
                    limitVal::pAdaptive::limitFlagLocal = true;
                }
                output[0] = 0.5*(UBCMod - UADMod);
            }
            else
            {
                output[0] = elementConsValOfOrder[1];
            }

            if (mathVar::orderElem > 1)
            {
                inputArgument[0] = UCD;
                inputArgument[1] = JPlus - elementConsValOfOrder[0];
                inputArgument[2] = elementConsValOfOrder[0] - JMinus;
                UCDMod = limiter::mathForLimiter::minmod(inputArgument);
                UCD_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);

                inputArgument[0] = -UAB;
                UABMod = -limiter::mathForLimiter::minmod(inputArgument);
                inputArgument[0] = UAB;
                UAB_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);
                if (fabs(UCD - UCD_check)>1e-10 || fabs(UAB - UAB_check)>1e-10)
                {
                    //std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
                    output[1] = 0.5*(UCDMod - UABMod);
                }
                else
                {
                    output[1] = elementConsValOfOrder[2];
                }
            }
            return output;
        }

        std::vector<double> pAdaptiveChildFunction_Tri(int element, int valType, double IPlus, double IMinus, double JMinus)
        {
            std::vector<double>elementConsValOfOrder(mathVar::orderElem + 1, 0.0),
                inputArgument(3, 0.0), output(mathVar::orderElem, 0.0);
            elementConsValOfOrder = auxUlti::getElementConserValuesOfOrder(element, valType);

            double UBC(math::pointValueNoLimiter(element, 1.0, 0.0, valType) - elementConsValOfOrder[0]),
                UAD(math::pointValueNoLimiter(element, -1.0, 0.0, valType) - elementConsValOfOrder[0]),
                UAB(math::pointValueNoLimiter(element, 0.0, -1.0, valType) - elementConsValOfOrder[0]),
                UBCMod(0.0), UADMod(0.0), UABMod(0.0),
                UBC_check(0.0), UAD_check(0.0), UAB_check(0.0),
                M(limiter::mathForLimiter::calM(element, valType)), Lxy(meshVar::localCellSize[element]);

            inputArgument[0] = -UAB;
            inputArgument[1] = elementConsValOfOrder[0] - JMinus;
            inputArgument[2] = elementConsValOfOrder[0] - JMinus;
            UABMod = -limiter::mathForLimiter::minmod(inputArgument);
            //inputArgument[0] = UAB;
            UAB_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);
            if (fabs(-UAB - UAB_check)>1e-10)
            {
                //std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
                output[0] = -UABMod;
            }
            else
            {
                output[0] = elementConsValOfOrder[1];
            }

            if (mathVar::orderElem > 1)
            {
                inputArgument[0] = UBC;
                inputArgument[1] = IPlus - elementConsValOfOrder[0];
                inputArgument[2] = elementConsValOfOrder[0] - IMinus;
                UBCMod = limiter::mathForLimiter::minmod(inputArgument);
                UBC_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1); //*fabs(elementConsValOfOrder[0])

                inputArgument[0] = -UAD;
                UADMod = -limiter::mathForLimiter::minmod(inputArgument);
                //inputArgument[0] = UAD;
                UAD_check = limiter::mathForLimiter::modifiedMinmod(inputArgument, M*Lxy*Lxy + 0.1);  //*fabs(elementConsValOfOrder[0])
                if ((fabs(UBC - UBC_check)>1e-10 || fabs(-UAD - UAD_check)>1e-10))
                {
                    //std::cout << "p-adaptive limiter is applied at cell " << element << " for variable type " << valType << std::endl;
                    if (limitVal::pAdaptive::limitFlagLocal == false)
                    {
                        limitVal::pAdaptive::limitFlagLocal = true;
                    }
                    output[1] = (UBCMod - UADMod) / 2.0;
                }
                else
                {
                    output[1] = elementConsValOfOrder[2];
                }
            }
            return output;
        }
    }
}
