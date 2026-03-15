#ifndef U_INTERIORSIDE_H
#define U_INTERIORSIDE_H
#include <fstream>
#include <vector>

namespace interiorSide {
    /**
     * @brief Function reads input of reflectRhoGrad condition.
     * @param bcGrp: group ID of boundary.
     * @param FileFlux: file flux to read.
     */
    void u_IO(int bcGrp, std::ifstream &FileFlux);

    /**
     * @brief Function corrects velocity following interiorSide condition (using zeroGradient_vector function).
     * @param varM: minus side velocity.
     * @param varP: plus side velocity.
     */
    void correctU(std::vector<double> &varM, std::vector<double> varP);

    /**
     * @brief Function corrects gradient of velocity following interiorSide condition.
     *
     * gradM = gradP
     *
     * @param gradM: minus side grad(U)
     * @param gradP: plus side grad(U)
     */
    void correctGradU(std::vector<double> &gradM, const std::vector<double> &gradP);
}

#endif // U_INTERIORSIDE_H
