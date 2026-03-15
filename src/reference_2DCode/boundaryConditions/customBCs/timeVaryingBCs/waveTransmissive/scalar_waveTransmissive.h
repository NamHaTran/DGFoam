#ifndef SCALAR_WAVETRANSIMISSIVE_H
#define SCALAR_WAVETRANSIMISSIVE_H
#include <fstream>
#include <vector>
namespace waveTransmissive {
    void solveScalarEqn(int edge, int edgeGrp, int nG, std::string var);

    void solveScalarEqn_implicit(int edge, int edgeGrp, int nG, std::string var);

    void calcDivPhi_DG(double *divPhi, int edge, int element, int nG, std::string var);

    void calcDivPhi_FDM(double *divPhi, int edge, int element, int nG, std::string var);

    //void correctPhi(int edge, int edgeGrp, int nG, std::vector<double> &varM, std::vector<double> varP);
}

#endif // SCALAR_WAVETRANSIMISSIVE_H
