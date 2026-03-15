#ifndef DGLIMITERLIB_H_INCLUDED
#define DGLIMITERLIB_H_INCLUDED
#include <tuple>
#include <vector>
namespace limiter
{
	//Function applies limiter to input element
    void limiter_1Step();

    void limitRho_PositivityPreserving();

	namespace Pp
	{
        void limiter();

		//Set initial value for theta1 and theta2
		void initialiseThetaVector();

        double calcTheta1(int element);

        double calcTheta2(int element);
	}

	namespace pAdaptive
	{
        void limiter();

        void limiter_1Elem(int element, int valType);

		std::vector<double> pAdaptiveChildFunction_Quad(int element, int valType, double IPlus, double IMinus, double JPlus, double JMinus);

		std::vector<double> pAdaptiveChildFunction_Tri(int element, int valType, double IPlus, double IMinus, double JMinus);
	}

    namespace massDiffusionLimiter
    {
        void limiter();

        bool checkTroubleCells();

        bool checkRunning();

        namespace mathForMassDiffLimiter {
            void calcVolumeGaussRho(int element);

            void calcRhoAtElementEdge(int element);

            void solveDivRho(int element);

            void calcFinvFvisAtInterface(int element);

            std::tuple<double, double, double, double> calcFinvFvisInVolume(int element, int na, int nb);

            void calcVolumeIntegralTerms(int element, std::vector<double>&VolIntTerm1);

            void calcSurfaceIntegralTerms(int element, std::vector<double>&SurfIntTerm1);

            void solveMassEquation(int RKOrder, int element);

            void updateVariables();

            double pointRho0(int element, double a, double b);

            double calcMinRhoResidual(int element);

            void limiter_1step(int RKOrder);
        }
    }

	namespace mathForLimiter
	{
        namespace Pp {
            double calcMinRho(int element);

            std::tuple<double, double> calcMinMeanRhoe(int element, double theta1);

            double calcRhoModified(int element, double a, double b, double theta1);

            double calRhoeFromConserVars(double rho, double rhou, double rhov, double rhoE);

            double calcTheta1Coeff(double minRho, double meanRho);

            double calcTheta2Coeff(double meanRhoe, double minRhoe, double meanRho);
        }

		//Function gets sign of input double number (return -1, 0, 1 if input number is negative, zero or positive, respectively)
		int getSignOfDouble(double input);

		//Minmod function
		double minmod(std::vector<double> vectorArguments);

		//Modified minmon function
		double modifiedMinmod(std::vector<double> vectorArguments, double condition);

		double calM(int element, int valType);

		void getNeighborElements();
	}
}

namespace troubledCellDetection
{
	bool checkTroubledCell(std::vector<double> InputVector, double condition);

    bool checkTroubleCell_massDiff(int elem);
}
#endif // DGLIMITERLIB_H_INCLUDED
