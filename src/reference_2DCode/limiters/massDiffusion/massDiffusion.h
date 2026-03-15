#ifndef MASSDIFFUSION_H
#define MASSDIFFUSION_H
#include <tuple>
#include <vector>

namespace limiter
{
    namespace massDiffusion
    {
        //Variables
        extern double DmCoeff, pseudoCo;

        extern bool trbCellAtMatchedBC;

        //! Max iteration of mass diffusion limiter
        extern int maxIter;

        //! Vector contents status (true/false) of all element which has at least 1 "matched BC type" edge.
        extern bool *markerOfTrbCellAtMatchedBC;
        extern bool *markerOfTrbCellAtMatchedBC_buffer;

        //! Counter of trouble edges which also be matched bc.
        extern int numOfTrbEdge;

        //Debug
        extern int RKOrder;


        void limiter();

        bool checkTroubleCells();

        std::tuple<bool, int> checkRunning();

        void updateVariables();

        void limitRho_PPLimitiedVer();

        namespace mathFuncs
        {
            void limiter_1step(int RKOrder);

            void CalcRHSTermOfDivRhoEqn(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy);

            void calcVolumeIntegralTermsOfDivRhoEqn(int element, std::vector<double> &rhoVolIntX, std::vector<double> &rhoVolIntY);

            void calcSurfaceIntegralTermsOfDivRhoEqn(int element, std::vector<double> &rhoSurfIntX, std::vector<double> &rhoSurfIntY);

            void getGaussVectorOfRho(int element, std::vector<std::vector<double>> &rhoFluxX, std::vector<std::vector<double>> &rhoFluxY);

            void calcVolumeGaussRho();

            void calcRhoAtElementSurf();

            void solveDivRho();

            void calcSurfaceConvDiffFlux();

            std::tuple<double, double, double, double> calcGaussConvDiffTerms(int element, int na, int nb);

            void calcVolumeIntegralConvDiffTerms(int element, std::vector<double>&VolIntTerm1);

            void calcSurfaceIntegralTerms(int element, std::vector<double>&SurfIntTerm1);

            void solveMassEquation(int RKOrder);

            void calcVolIntegralTerm(int element,std::vector<double>&VolIntTerm1);

            double calcMinRhoResidual(int element);

            double pointRho0(int element, double a, double b);

            double calcMinUWOLimiter(int element, int type);

            double calcMaxAbsGradRhoOnEdges(int element);

            double calcMeanAbsGradRho(int element);
        }
    }

    namespace IOMassDiff
    {
        void readSetting();
    }
}

#endif // MASSDIFFUSION_H
