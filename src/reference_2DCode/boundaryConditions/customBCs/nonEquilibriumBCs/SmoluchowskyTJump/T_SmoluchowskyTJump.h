#ifndef T_SMOLUCHOWSKYTJUMP_H
#define T_SMOLUCHOWSKYTJUMP_H
#include <fstream>
#include <vector>

namespace SmoluchowskyTJump {
    /**
     * @brief Function reads input of Smolukowsky temperature jump boundary condition.
     *
     * @param bcGrp: group ID of boundary.
     * @param FileFlux: file flux to read.
     */
    void T_IO(int bcGrp, std::ifstream &FileFlux);

    void correctT(int edge, int edgeGrp, int nG, double &varM, double varP);

    //void correctGradT(std::vector<double> &gradM, const std::vector<double> &gradP);

    /**
     * @brief Function calculates jumped temperature at wall following Smolukowsky's condition.
     *
     * This function uses explicit method, in which gradient of T is calculated explicitly from cell's data using function math::pointAuxValue.
     *
     * All data are calculated at a normal projection of cell's centroid to considering edge, so this function runs only 1 time for 1 edge.
     *
     * Smoluchowsky temperature jump equation (refer in "Langmuir–Maxwell and Langmuir–Smoluchowski boundary conditions for thermal gas flow simulations in hypersonic aerodynamics", Nam T.P. Le et al, ):
     *
     * \f$ T + \left( \frac{ 2 - \sigma_T }{ \sigma_T } \right) \frac{ 2 \gamma }{ \left( \gamma + 1 \right) Pr } \lambda \nabla_n T = T_w \f$.
     *
     * with
     *
     * - \f$ \lambda = \frac{ \mu }{ \rho } \sqrt{ \frac{ \pi }{ R T } } \f$ is mean free path.
     *
     * @param edge: edge ID needed to apply Maxwell's condition.
     * @param edgeGrp: group ID of edge.
     * @param nG: ID of Gauss point.
     */
    void calcTJump_DGTypeExplicit(int edge, int edgeGrp, int nG);

    /**
     * @brief Function calculates jumped temperature at wall following Smolukowsky's condition.
     *
     * This function uses implicit method bases on FDM, in which gradient of T are approximated by using forward difference.
     *
     * All data are calculated at a normal projection of cell's centroid to considering edge, so this function runs only 1 time for 1 edge.
     *
     * Smoluchowsky temperature jump equation (refer in "Langmuir–Maxwell and Langmuir–Smoluchowski boundary conditions for thermal gas flow simulations in hypersonic aerodynamics", Nam T.P. Le et al, ):
     *
     * \f$ T + \left( \frac{ 2 - \sigma_T }{ \sigma_T } \right) \frac{ 2 \gamma }{ \left( \gamma + 1 \right) Pr } \lambda \nabla_n T = T_w \f$.
     *
     * with
     *
     * - \f$ \lambda = \frac{ \mu }{ \rho } \sqrt{ \frac{ \pi }{ R T } } \f$ is mean free path.
     *
     * Terms \f$ \nabla T \f$ is approximated using forward difference.
     *
     * - \f$ \nabla T = \frac{T_{jump} - T_C}{d} \f$ with \f$ T_{jump} \f$ and \f$ T_C \f$ are jump temperature at wall (unknow) and temperature at cell's centroid, \f$ d \f$ is normal distance from cell's centroid to wall.
     *
     * Note that in the function, viscosity coefficient is calculated using Sutherland's formula: \f$ \mu = A_s \frac{ T^{1.5} }{ T + T_s } \f$ .
     *
     * This equation and approximated \f$ \nabla T \f$ are substitued to Smoluchowsky's equation to get final equation of only 1 unknow is \f$ T_{jump} \f$ (\f$ T_{jump} \f$ is solved completely implicitly).
     *
     * Because of that, this function is not suitable for other viscosity models (such as powerLaw_VHS).
     *
     * To improve this, another function can be created in which \f$ \mu \f$ is calculated explicitly using inputted viscosity model.
     *
     * @param edge: edge ID needed to apply Smoluchowsky's condition.
     * @param edgeGrp: group ID of edge.
     * @param nG: ID of Gauss point.
     */
    void calcTJump_FDMTypeImplicit(int edge, int edgeGrp, int nG);

    void calcTJump_FDMTypeSemiImplicit(int edge, int edgeGrp, int nG);
}
#endif // T_SMOLUCHOWSKYTJUMP_H
