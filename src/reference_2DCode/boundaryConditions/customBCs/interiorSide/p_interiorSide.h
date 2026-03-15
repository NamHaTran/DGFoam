#ifndef P_INTERIORSIDE_H
#define P_INTERIORSIDE_H
#include <fstream>
#include <vector>
namespace interiorSide {
    /**
     * @brief Function reads input of reflectRhoGrad condition.
     * @param bcGrp: group ID of boundary.
     * @param FileFlux: file flux to read.
     */
    void p_IO(int bcGrp, std::ifstream &FileFlux);

    /**
     * @brief Function corrects pressure following interiorSide condition (using zeroGradient_scalar function).
     * @param pM: minus side pressure
     * @param pP: plus side pressure
     */
    void correctP(double &pM, double pP);

    /**
     * @brief Function corrects gradient of pressure following interiorSide condition.
     *
     * gradM = gradP
     *
     * @param gradM: minus side grad(p)
     * @param gradP: plus side grad(p)
     */
    void correctGradP(std::vector<double> &gradM, const std::vector<double> &gradP);
}

#endif // P_INTERIORSIDE_H
