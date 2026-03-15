#ifndef DGPROCLIB_H_INCLUDED
#define DGPROCLIB_H_INCLUDED
#include <vector>

namespace meshParam
{
    /**
     * @brief Function generates base vector of coordinates of Gauss points
     */
    void genBasedGaussPtsVectors();

    /**
     * @brief Function calculates coordinates and weights of Gauss Points
     */
	void GaussParam();

    /**
     * @brief Function calculates Jacobian of 2D element at all volume Gauss points
     */
	void JacobianParam();

    /**
     * @brief Function calculate basis function at all volume Gauss points
     */
	void basisFcParam();

    /**
     * @brief Function saves coordinates derivatives at all volume Gauss points to array
     */
	void derivCoordinates();

	/*Function calculates centroid and cell size of elements*/
	void calcCellMetrics();

    void calcEdgeLength();

	void calcStiffMatrixCoeffs();

    void calcDistanceFromCenterToBCEdge();

    void findNormProjectionOfCenterToBCEdge();
}

namespace process
{
    /*Ham chuan bi case truoc khi chay*/
    void prepareCaseForRunning();

    /*Ham lam cac thu tuc truoc khi ket thuc 1 iteration*/
    void finalizeIteration();

	/*Function sets initial value to all elements*/
    void setIniVolumeValues();

    /*Function sets initial value to all Gauss points at BC*/
    void setIniSurfaceBCValues();

	/*Function distributes initial value to each order of accuracy of element*/
	std::vector<double> calcIniValues(double iniVal, int element);

	/*Function computes RHS of initial condition equation*/
	std::vector<double> calcIniValuesRHS(int element, double iniVal);

    /**
     * @brief Function calculates value of conservative variables on cell's surfaces
     */
    void calcValuesAtInterface();

    /**
     * @brief Function calculates conservative variables at all volume Gauss points
     */
	void calcVolumeGaussValues();

    void calcVolumeGaussRho();

    void calcRhoAtInterface();

    /**
     * @brief Function calculates temperature at all Gauss points
     */
    void calcTGauss();

	namespace auxEq
	{
        /*Function solve div(rho) in case of mass diffusion*/
        void solveDivRho();

        /*Function solves auxilary equation at all elements for auxilary variables*/
        void solveAuxEquation();

        /*Ham nay su dung khi flow inviscid, de tinh toan cac gia tri U tai surfaceFields*/
        void updateSurfaceFieldsAtBC();

        namespace BR1 {
        void calcValuesAtInterface();

        std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType);

        /*Function calculates right hand side term of auxilary equations at all order*/
        void CalcRHSTerm(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy, std::vector<double> &rhouRHSOx, std::vector<double> &rhouRHSOy, std::vector<double> &rhovRHSOx, std::vector<double> &rhovRHSOy, std::vector<double> &rhoERHSOx, std::vector<double> &rhoERHSOy);

        /*Function calculates volume integral terms of element*/
        void calcVolumeIntegralTerms(int element, std::vector<double> &rhoVolIntX, std::vector<double> &rhouVolIntX, std::vector<double> &rhovVolIntX, std::vector<double> &rhoEVolIntX, std::vector<double> &rhoVolIntY, std::vector<double> &rhouVolIntY, std::vector<double> &rhovVolIntY, std::vector<double> &rhoEVolIntY);

        /*Function calculates surface integral terms of element*/
        void calcSurfaceIntegralTerms(int element, std::vector<double> &rhoSurfIntX, std::vector<double> &rhouSurfIntX, std::vector<double> &rhovSurfIntX, std::vector<double> &rhoESurfIntX, std::vector<double> &rhoSurfIntY, std::vector<double> &rhouSurfIntY, std::vector<double> &rhovSurfIntY, std::vector<double> &rhoESurfIntY);

        /*Function returns matrix content Gauss values of <valType> conservative variable at all Gauss points in element volume,
        use only for auxilary equation*/
        void getGaussMatrixOfConserVar(std::vector<std::vector<double>> &UGauss, int element, int valType);

        /*Function returns vector of flux of all conservative variables at all faces of element*/
        void getGaussVectorOfConserVar(int element, std::vector<std::vector<double>> &rhoFluxX, std::vector<std::vector<double>> &rhouFluxX, std::vector<std::vector<double>> &rhovFluxX, std::vector<std::vector<double>> &rhoEFluxX, std::vector<std::vector<double>> &rhoFluxY, std::vector<std::vector<double>> &rhouFluxY, std::vector<std::vector<double>> &rhovFluxY, std::vector<std::vector<double>> &rhoEFluxY);

        /*Function returns vector content Gauss values of all conservative variable at all Gauss points on the edge,
        use only for auxilary equation*/
        std::vector<std::vector<double>> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nx, double ny);
        }

        namespace BR2 {
            void calcVolDivU(int element);

            void calcSurDivU(int element);
        }
	}

	namespace NSFEq
	{
        std::tuple<double, double> getInterfacesFluxes(int edge, int element, int nG, int mod, int valType);

        void solveNSFEquation_Euler();

		/*Function solves NSEF equation at all elements for conservative variables*/
        void solveNSFEquation_TVDRK(int RKOrder);

		/*Function calculates right hand side terms of all conservative variables at ONLY one order*/
		void CalcRHSTerm(int element, std::vector<double> &term1RHS, std::vector<double> &term2RHS, std::vector<double> &term3RHS, std::vector<double> &term4RHS);

		/*Function calculates Inviscid terms at Gauss point (a, b) and returns 2D matrix
		InviscidTerm 2D array has 4 row 2 column:
		- column 1: Ox direction
		- column 2: Oy direction*/
		std::vector<std::vector<double>> calcGaussInviscidTerm(int element, int na, int nb);

		/*Function calculates Viscous terms at Gauss point (a, b) and returns 2D matrix
		ViscousTerm 2D array has 4 row 2 column:
		- column 1: Ox direction
		- column 2: Oy direction*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, int na, int nb);

		/*Function calculates volume integral terms in NSF equation at ONLY ONE ORDER*/
		void calcVolumeIntegralTerms(int element, std::vector<double> &VolIntTerm1, std::vector<double> &VolIntTerm2, std::vector<double> &VolIntTerm3, std::vector<double> &VolIntTerm4);

		/*Function calculates surface integral terms in NSF equation at ONLY ONE ORDER*/
		void calcSurfaceIntegralTerms(int element, std::vector<double> &SurfIntTerm1, std::vector<double> &SurfIntTerm2, std::vector<double> &SurfIntTerm3, std::vector<double> &SurfIntTerm4);

		/*Function calculates flux at nGauss point of all conservative variables at internal egde*/
        std::vector<double> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, int nGauss);

		/*Function applies input ddtSchemme to solve time marching*/
		std::vector<double> solveTimeMarching(int element, std::vector<double> &ddtArr, std::vector<double> &UnArr, int RKOrder, int varType);

        void calcLocalInviscidFlux(std::vector<double> &UG, std::vector<double> &FLocalX, std::vector<double> &FLocalY, std::vector<double> &UL, std::vector<double> &n, double T);

        void calcLocalViscousFlux(std::vector<double> &UG, std::vector<double> &dUX, std::vector<double> &dUY, std::vector<double> &FLocalX, std::vector<double> &FLocalY, std::vector<double> &n);

        void calcSurfaceFlux();
	}

	namespace timeDiscretization
	{
        void solveTimeEq();

        /**
         * @brief Function solves 1 time steps of DG workflow by using Euler time discretization
         */
        void Euler();

		void calcGlobalTimeStep();

		/*Function computes local time step*/
		double localTimeStep(int element);

		/*Function computes time derivative at centroid of standard cell*/
		double localErrorEstimate(int element, std::vector<double> &ddtArr);

		/*Function computes global errors*/
		void globalErrorEstimate();

		void TVDRK_1step(int RKOrder);

		void TVDRK3();

        namespace parallel {
        void minTimeStep();
        }
	}

	/*Function calculates volume integral terms of auxilary equations of ONLY one order
	User's guide:
	elem: element index
	Ui: array of Gauss values
	direction: index of direction
	1: x direction
	2: y direction*/
	double volumeInte(int elem, std::vector< std::vector<double> > &Ui, int order, int direction);

	/*Function calculates surface integral at 1 surface of auxilary equations of ONLY one order
	User's guide:
	direction: index of direction
	1: x direction
	2: y direction*/
	double surfaceInte(int elem, int edge, std::vector<double> &FluxVector, int order);

	/*Function calculates Auxilary Stiff matrix*/
	std::vector<std::vector<double>> calculateStiffMatrix(int element);

    /**
     * @brief Function calculates an element of stiff matrix (volume integral)
     * @param element: element ID
     * @param order1: order 1
     * @param order2: order 2
     * @return value of the element of stiff matrix
     */
	double calculateStiffMatrixElement(int element, int order1, int order2);

	/*Function return false if running contion is wrong*/
	bool checkRunningCond();

	std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType);
}

#endif // DGPROCLIB_H_INCLUDED
