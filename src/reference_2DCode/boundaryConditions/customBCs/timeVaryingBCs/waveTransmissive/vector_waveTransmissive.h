#ifndef VECTOR_WAVETRANSIMISSIVE_H
#define VECTOR_WAVETRANSIMISSIVE_H
#include <fstream>
#include <vector>
namespace waveTransmissive {
    void solveVectorEqn(int edge, int edgeGrp, int nG);

    void calcDivU_DG(double *divU, int edge, int element, int nG);

    void calcDivU_FDM(double *divU, int edge, int element, int nG);

    //void correctU(int edge, int edgeGrp, int nG, std::vector<double> &varM, std::vector<double> varP);
}

#endif // VECTOR_WAVETRANSIMISSIVE_H
