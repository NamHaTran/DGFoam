#ifndef T_INTERIORSIDE_H
#define T_INTERIORSIDE_H
#include <fstream>
#include <vector>

namespace interiorSide {
    /**
     * @brief Function reads input of reflectRhoGrad condition.
     * @param bcGrp: group ID of boundary.
     * @param FileFlux: file flux to read.
     */
    void T_IO(int bcGrp, std::ifstream &FileFlux);

    /**
     * @brief Function corrects temperature following interiorSide condition (using zeroGradient_scalar function).
     * @param varM: minus side temperature.
     * @param varP: plus side temperature.
     */
    void correctT(double &varM, double varP);

    /**
     * @brief Function corrects gradient of temperature following interiorSide condition.
     *
     * gradM = gradP
     *
     * @param gradM: minus side grad(T)
     * @param gradP: plus side grad(T)
     */
    void correctGradT(std::vector<double> &gradM, const std::vector<double> &gradP);
}
#endif // T_INTERIORSIDE_H
