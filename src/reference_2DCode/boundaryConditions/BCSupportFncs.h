#ifndef BCSUPPORTFNCS_H
#define BCSUPPORTFNCS_H
#include <vector>
namespace BCSupportFncs
{
    void calcBCValuesAtUnfixedValueBCs();

    /**
     * @brief Function checks inflow/outflow.
     *
     * Function returns true if at considering face (edge), flow is going into computational domain, returns false if flow is going out of computational domain.
     *
     * @param u: Ox velocity component.
     * @param v: Oy velocity component.
     * @param nx: Ox component of normal unit vector.
     * @param ny: Oy component of normal unit vector.
     * @return
     */
    bool checkInflow(double u, double v, double nx, double ny);

    /**
     * @brief Function decomposes conservative variables (U) to primary variables.
     * @param priVars: vector of primary variables (returned value).
     * @param U: vector of conservative variables.
     * @param calcT: flag of calculating T from conservative variables (detail is in function).
     */
    void decompseU(std::vector<double> &priVars, const std::vector<double> &U, bool calcT);

    /**
     * @brief Function decomposes gradient of conservative variables (U) to gradient of primary variables.
     * @param priVars: vector of gradient of primary variables (returned value).
     * @param U: vector of conservative variables.
     * @param dU: vector of gradient of conservative variables.
     * @param T: Temperature.
     */
    void decompsedU(std::vector<double> &priVars, const std::vector<double> &U, const std::vector<double> &dU, double T);

    /**
     * @brief Function correct primary variables at minus side by using variables at plus side and boundary conditions.
     * @param edge: edge ID.
     * @param edgeGrp: group ID of edge.
     * @param nG: ID of Gauss point.
     * @param priVarsM: vector of primary variables at minus side (returned value).
     * @param priVarsP: vector of primary variables at plus side.
     * @param priVarsMean: vector of mean values of primary variables at plus side.
     * @param n: normal unit vector.
     * @param inflow: flag of inflow.
     */
    void correctPriVars(int edge, int edgeGrp, int nG, std::vector<double> &priVarsM, const std::vector<double> &priVarsP, const std::vector<double> &priVarsMean, const std::vector<double> &n, bool inflow);

    /**
     * @brief Function correct gradient of primary variables at minus side by using variables at plus side and boundary conditions.
     * @param edge: edge ID.
     * @param edgeGrp: group ID of edge.
     * @param nG: ID of Gauss point.
     * @param dpriVarsXM: vector of gradient of primary variables at minus side on Ox direction (returned value).
     * @param dpriVarsYM: vector of gradient of primary variables at minus side on Oy direction (returned value).
     * @param dpriVarsXP: vector of gradient of primary variables at plus side on Ox direction.
     * @param dpriVarsYP: vector of gradient of primary variables at plus side on Oy direction.
     * @param UP: vector of conservative variables at plus side.
     * @param UM: vector of conservative variables at minus side.
     * @param TP: temperature at plus side.
     * @param TM: temperature at minus side.
     * @param n: normal unit vector.
     * @param inflow: flag of inflow.
     */
    void correctPriVarsGrad(int edge, int edgeGrp, int nG, std::vector<double> &dpriVarsXM, std::vector<double> &dpriVarsYM, const std::vector<double> &dpriVarsXP, const std::vector<double> &dpriVarsYP, const std::vector<double> &UP, const std::vector<double> &UM,  double TP, double TM, const std::vector<double> &n, bool inflow);

    /**
     * @brief Function correct density at minus side by using density at plus side and boundary conditions.
     * @param priVarsM: density at minus side (returned value).
     * @param edgeGrp: group ID of edge.
     * @param priVarsP: density at plus side.
     */
    void correctDensity(std::vector<double> &priVarsM, int edgeGrp, const std::vector<double> &priVarsP);

    /**
     * @brief Function correct gradient of density at minus side by using density at plus side and boundary conditions.
     * @param edgeGrp: group ID of edge.
     * @param dpMX: gradient of pressure at minus side on Ox direction.
     * @param dpMY: gradient of pressure at minus side on OY direction.
     * @param dTMX: gradient of temperature at minus side on Ox direction.
     * @param dTMY: gradient of temperature at minus side on Oy direction.
     * @param dRhoPX: gradient of density at plus side on Ox direction (returned value).
     * @param dRhoPY: gradient of density at plus side on Oy direction (returned value).
     * @param UM: vector of conservative variables at minus side.
     * @param TM: temperature at minus side.
     * @param n: normal unit vector.
     * @return 2 components of div(density) on minus side.
     */
    std::tuple<double, double> correctDensityGrad(int edgeGrp, double dpMX, double dpMY, double dTMX, double dTMY, double dRhoPX, double dRhoPY, const std::vector<double> &UM, double TM, const std::vector<double> &n);

    std::tuple <double, double> correctVelocity(int edge, int edgeGrp, int nG, double uP, double uMean, double vP, double vMean, const std::vector<double> &n, bool inflow);

    std::tuple <double, double> correctVelocityGrad(int edgeGrp, int nG, double duXP, double duYP, bool inflow, const std::vector<double> &n);

    double correctTemperature(int edge, int edgeGrp, int nG, double TP, double TMean, bool inflow);

    std::tuple<double, double> correctTemperatureGrad(int edgeGrp, int nG, double dTXP, double dTYP, bool inflow, const std::vector<double> &n);

    double correctPressure(int edge, int edgeGrp, int nG, double pP, double pMean, bool inflow);

    std::tuple<double, double> correctPressureGrad(int edgeGrp, double dpXP, double dpYP, bool inflow, const std::vector<double> &n);

    void reconstructConvectiveU(std::vector<double> &U, const std::vector<double> &priVars);

    void reconstructdU(std::vector<double> &dU, const std::vector<double> &priVars, const std::vector<double> &U, double T);

    namespace auxilaryBCs {
        void getUPlus(int edge, int nG, std::vector<double> &UPlus);

        void calcUMean(int element, std::vector<double> &UMean);

        /*Ham apply dieu kien back flow tai vi tri bi back (reversed) flow*/
        void calcUReversedFlow(int edge, int nG, std::vector<double> &U);
    }

    namespace NSFEqBCs {
        void findMaxLxFConstantOnEdge(int element, int edge);

        void getdUPlus(int edge, int nG, std::vector<double> &dUXPlus, std::vector<double> &dUYPlus);

        void calcdUMean(int edge, int element, std::vector<double> &dUXMean, std::vector<double> &dUYMean);

        /*Function calculates fluxes of all advective and diffusive terms of NSF equation from conservative variables*/
        void NSFEqFluxes(int edgeId, std::vector<double> &Fluxes, double TPlus, double TMinus, std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector);
    }

    namespace parallel {
    void getUMinus_convective(int edge, int nG, std::vector<double> &UMinus);

    void getdUMinus(int edge, int nG, std::vector<double> &dUMinus, int dir);
    }
}
#endif // BCSUPPORTFNCS_H
