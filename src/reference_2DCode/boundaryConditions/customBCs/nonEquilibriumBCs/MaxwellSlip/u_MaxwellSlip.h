#ifndef U_MAXWELLSLIP_H
#define U_MAXWELLSLIP_H
#include <fstream>
#include <vector>

namespace MaxwellSlip {

    /**
     * @brief Function reads input of Maxwell slip boundary condition.
     *
     * @param bcGrp: group ID of boundary.
     * @param FileFlux: file flux to read.
     */
    void u_IO(int bcGrp, std::ifstream &FileFlux);

    void correctU(int edge, int edgeGrp, int nG, std::vector<double> &varM, std::vector<double> varP);

    //void correctGradU(std::vector<double> &gradM, const std::vector<double> &gradP);

    /**
     * @brief Function calculates slip velocity at wall following Maxwell's condition.
     *
     * This function uses explicit method, in which gradients of U(u,v) and T are calculated explicitly from cell's data using function math::pointAuxValue.
     *
     * All data are calculated at a normal projection of cell's centroid to considering edge, so this function runs only 1 time for 1 edge.
     *
     * Maxwell slip equation (inclued curvature effect and thermal creep, refer in "Langmuir–Maxwell and Langmuir–Smoluchowski boundary conditions for thermal gas flow simulations in hypersonic aerodynamics", Nam T.P. Le et al, ):
     *
     * \f$ \textbf{u} + \left( \frac{ 2 - \sigma_u }{ \sigma_u } \right) \lambda \nabla_n \left( \textbf{S} \cdot \textbf{u} \right) = \textbf{u}_w - \left( \frac{2 - \sigma_u}{\sigma_u} \right) \frac{\lambda}{\mu} \textbf{S} \cdot \left( \textbf{n} \cdot \Pi_{mc} \right) - \frac{3}{4} \frac{ \mu}{ \rho} \frac{ \textbf{S} \cdot \nabla T }{ T }\f$
     *
     * with
     *
     * - Term \f$ \Pi_{mc} = \mu \left( \nabla \textbf{u}^{T} \right) - \frac{2}{3} \textbf{I} tr \left( \nabla \textbf{u} \right)\f$
     *
     *
     * @param edge: edge ID needed to apply Maxwell's condition.
     * @param edgeGrp: group ID of edge.
     * @param nG: ID of Gauss point.
     */
    void calcUSlip_DGTypeExplicit(int edge, int edgeGrp, int nG);

    /**
     * @brief Function calculates slip velocity at wall following Maxwell's condition.
     *
     * This function uses implicit method bases on FDM, in which gradients of U(u,v) and T are approximated by using forward difference.
     *
     * All data are calculated at a normal projection of cell's centroid to considering edge, so this function runs only 1 time for 1 edge.
     *
     * Maxwell slip equation (inclued curvature effect and thermal creep, refer in "Langmuir–Maxwell and Langmuir–Smoluchowski boundary conditions for thermal gas flow simulations in hypersonic aerodynamics", Nam T.P. Le et al, ):
     *
     * \f$ \textbf{u} + \left( \frac{ 2 - \sigma_u }{ \sigma_u } \right) \lambda \nabla_n \left( \textbf{S} \cdot \textbf{u} \right) = \textbf{u}_w - \left( \frac{2 - \sigma_u}{\sigma_u} \right) \frac{\lambda}{\mu} \textbf{S} \cdot \left( \textbf{n} \cdot \Pi_{mc} \right) - \frac{3}{4} \frac{ \mu}{ \rho} \frac{ \textbf{S} \cdot \nabla T }{ T }\f$.
     *
     * with
     *
     * - Term \f$ \Pi_{mc} = \mu \left( \nabla \textbf{u}^{T} \right) - \frac{2}{3} \textbf{I} tr \left( \nabla \textbf{u} \right)\f$
     *
     * Terms \f$ \nabla T \f$ and \f$ \nabla (u, v) \f$ are approximated using forward difference.
     *
     * - \f$ \nabla T = \frac{T_{jump} - T_C}{d} \f$ with \f$ T_{jump} \f$ and \f$ T_C \f$ are jump temperature at wall and temperature at cell's centroid, d is normal distance from cell's centroid to wall.
     *
     * - \f$ \nabla (u, v) = \frac{(u, v)_{slip} - (u, v)_C}{d} \f$ with \f$ (u, v)_{slip} \f$ are slip velocity components (unknows) and \f$ (u, v)_C \f$ velocity at cell's centroid.
     *
     * Substitue approximated \f$ \nabla (u, v) \f$ to Maxwell's slip equation and solve equation to get \f$ (u, v)_{slip} \f$.
     *
     * @param edge: edge ID needed to apply Maxwell's condition.
     * @param edgeGrp: group ID of edge.
     * @param nG: ID of Gauss point.
     */
    void calcUSlip_FDMTypeImplicit(int edge, int edgeGrp, int nG);

    double calcGradientUsingFDM(double phiC, double phiBC, double delta, double idotn);
}

#endif // U_MAXWELLSLIP_H
