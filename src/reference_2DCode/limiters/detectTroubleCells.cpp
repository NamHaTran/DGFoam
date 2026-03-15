#include "detectTroubleCell.h"
#include <vector>
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include "dynamicVarDeclaration.h"
#include <math.h>
#include "mathFunctions.h"

namespace troubledCellDetector
{
    bool minmodDetector(std::vector<double> InputVector, double condition)
    {
        bool needLimit(true);
        if (fabs(limiter::mathForLimiter::modifiedMinmod(InputVector,condition) - InputVector[0]) > 1e-10)
        {
            needLimit = false;
        }
        else
        {
            needLimit = true;
        }
        return needLimit;
    }

    /**
     * @brief Function detects discontinuity of density solution bases on Persson & Paraire's shock detection algorithm.
     *
     * Persson & Paraire's shock detection algorithm is based on "the rate of decay of the expansion coefficients of the solution when this is expressed in a hierarchical orthonormal basis".
     *
     * Reference: P. Persson and J. Peraire, “Sub-cell shock capturing for discontinuous Galerkin methods,” AIAA Paper 2006-0112, Jan. 2006.
     *
     * @param elem: element Id
     * @return
     */
    bool PerssonPeraireDetector(int elem)
    {
        //Dung discontinuity dectector (Persson & Peraire)
        bool troubleCell(false);
        std::vector<std::vector<double>>
            SnumGs(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0)),
            SdenGs(mathVar::nGauss2D + 1, std::vector<double>(mathVar::nGauss2D + 1, 0.0));
        double a(0.0), b(0.0), S(0.0);

        for (int na = 0; na <= mathVar::nGauss2D; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
            {
                std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
                math::basisFc(a, b, auxUlti::checkType(elem));  //mathVar::B is changed after this command line is excuted

                SnumGs[na][nb] = pow(mathVar::B[mathVar::orderElem]*rho[elem][mathVar::orderElem],2);
                SdenGs[na][nb]=0.0;
                for (int iorder=0; iorder<mathVar::orderElem; iorder++)
                {
                    SdenGs[na][nb]+=pow(mathVar::B[iorder]*rho[elem][iorder],2);
                }
                SdenGs[na][nb]+=SnumGs[na][nb];
            }
        }

        S=math::volumeInte(SnumGs,elem)/math::volumeInte(SdenGs,elem);
        //Check condition
        /* Condition is refered from Persson's paper: "we expect that the value of Se will scale like ∼ 1/p^4".
        */
        if (S>1/pow(mathVar::orderElem+1,6))
        {
            troubleCell=true;
        }
        return troubleCell;
    }
}
