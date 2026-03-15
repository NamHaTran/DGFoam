#ifndef DURSTMODEL_H
#define DURSTMODEL_H
#include <tuple>
#include <vector>

namespace extNSF_Durst {
    /*Variables*/
    extern bool
    enable,
    diffusionAtWall,
    needToRemoveDiffTerm,

    massDiffModel_ChapmanEnskog,
    massDiffModel_constant,

    dropNormSelfDiffTerm,

    isSymmetry;

    extern double Dm,
    blending,
    realDm;

    //mass diffusion velocity at wall
    extern double *uMassWall;
    extern double *vMassWall;

    void resizeArrays();

    void applyBlendingFactorToDm();

    void correctViscousTerms(std::vector<std::vector<double>> &diffTerms, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy);

    std::tuple<double, double, double, double> calcSelfDiffusionTerms(std::vector< std::vector<double> > &selfDiffusionTensor, double uVal, double vVal, int dir);

    std::vector<std::vector<double>> calcSelfDiffusionTensor(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy);

    double calcDiffusionStressComponent(int index, double fstTerm, double sndTerm);

    void correctEnergyEqnVolIntTerm(int element, std::vector<double> &VolIntTerm4);

    double calcSelfDiffFlux(double rho, double T, double gradRho, double gradT);

    double calcDiffVelocity(std::vector<double> &U, std::vector<double> &dU);
}

#endif // DURSTMODEL_H
