#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
#include <vector>

namespace bcForExtNSF_Durst {
extern double
//! Variable for correcting boundary condition in case of Durst model. tempdTx is temporary dTx.
tempdTx,
//! Variable for correcting boundary condition in case of Durst model. tempdTy is temporary dTy.
tempdTy;

//void dropNormDiffVel(int edge, int edgeGrp, double &dRhoXM, double &dRhoYM, double dRhoXP, double dRhoYP, double rhoBC, double dTBCx, double dTBCy, const std::vector<double> &n);

void checkConditionToAddDiffTerms(int edge);

void checkConditionToDropNormSelfDiffTerm(int edge);

void resetNeedToRemoveDiffTermFlag();

void resetNeedToDropNormSelfDiffTermFlag();

void correctViscousTerms(std::vector<std::vector<double>> &diffTerms, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy, std::vector<double> &n);

std::vector<std::vector<double>> calcSelfDiffusionTensor(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy, std::vector<double> &n);

void removeNormTermOfVector(double &dx, double &dy, const std::vector<double> &n, const double vectorDirX, const double vectorDirY, bool removeOnlyWhenVectorDirIsSameDirectionWithn);
}
#endif // BOUNDARYCONDITIONS_H
