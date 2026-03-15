#include "mathFunctions.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include "DGAuxUltilitiesLib.h"
#include "dynamicVarDeclaration.h"

namespace limiter
{
    //MATHEMATIC FUNCTIONS FOR LIMITER
    namespace mathForLimiter
    {
        double minmod(std::vector<double> vectorArguments)
        {
            int numberArgument(static_cast<int>(vectorArguments.size()));
            std::vector<double>absVectorArguments(numberArgument, 0.0);
            double minmodVal(0.0);
            bool isSameSign(true);
            for (int i = 0; i < numberArgument-1; i++)
            {
                if (limiter::mathForLimiter::getSignOfDouble(vectorArguments[i]) != limiter::mathForLimiter::getSignOfDouble(vectorArguments[i+1]))
                {
                    isSameSign = false;
                    break;
                }
                absVectorArguments[i] = fabs(vectorArguments[i]);
            }
            absVectorArguments[numberArgument - 1] = fabs(vectorArguments[numberArgument - 1]);

            if (isSameSign)
            {
                minmodVal = (*std::min_element(absVectorArguments.begin(), absVectorArguments.end())) * limiter::mathForLimiter::getSignOfDouble(vectorArguments[0]);
            }
            else
            {
                minmodVal = 0.0;
            }
            return minmodVal;
        }

        double modifiedMinmod(std::vector<double> vectorArguments, double condition)
        {
            double minmodVal(0.0);
            if (fabs(vectorArguments[0])<=condition)
            {
                minmodVal = vectorArguments[0];
            }
            else
            {
                minmodVal = limiter::mathForLimiter::minmod(vectorArguments);
            }
            return minmodVal;
        }

        int getSignOfDouble(double input)
        {
            int sign(0);
            if (input < 0)
            {
                sign = -1;
            }
            else if (input == 0.0)
            {
                sign = 0;
            }
            else
            {
                sign = 1;
            }
            return sign;
        }

        double calM(int element, int valType)
        {
            int neighborElemId(0);
            double M(0.0);
            for (int i = 0; i < auxUlti::checkType(element); i++)
            {
                neighborElemId = meshVar::esuel[element][i];
                if (neighborElemId<0)
                {
                    neighborElemId = element;
                }
                switch (valType)
                {
                case 1:
                {
                    M += rho[neighborElemId][0] - rho[element][0];
                }
                break;
                case 2:
                {
                    M += rhou[neighborElemId][0] - rhou[element][0];
                }
                break;
                case 3:
                {
                    M += rhov[neighborElemId][0] - rhov[element][0];
                }
                break;
                case 4:
                {
                    M += rhoE[neighborElemId][0] - rhoE[element][0];
                }
                break;
                default:
                    break;
                }
            }
            M = fabs(M);
            return M;
        }

        void getNeighborElements()
        {
            for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
            {
                for (int i = 0; i < auxUlti::checkType(nelem); i++)
                {
                    meshVar::neighboringElements[nelem][i] = auxUlti::getNeighborElement(nelem, auxUlti::getEdgeHasInputOrderOfElement(nelem, i));
                }
            }
        }
    }
}
