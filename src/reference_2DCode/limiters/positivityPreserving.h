#ifndef POSITIVITYPRESERVING_H
#define POSITIVITYPRESERVING_H

#include <tuple>

namespace limiter
{
    namespace positivityPreserving
    {
        void limiter();

        void initialiseThetaVector();

        double calcTheta2(int element);

        double calcTheta1(int element);

        namespace mathFuncs
        {
            double calcMinRho(int element);

            //Function calculates modified value of Rho at abitrary point (for calculating theta2)
            double calcRhoModified(int element, double a, double b, double theta1);

            std::tuple<double, double> calcMinMeanRhoe(int element, double theta1);

            double calRhoeFromConserVars(double rho, double rhou, double rhov, double rhoE);

            double calcTheta1Coeff(double minRho, double meanRho);

            double calcTheta2Coeff(double meanRhoe, double minRhoe, double meanRho);
        }
    }
}
#endif // POSITIVITYPRESERVING_H
