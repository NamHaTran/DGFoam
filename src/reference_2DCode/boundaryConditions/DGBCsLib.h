/* Huong dan add dieu kien bien moi vao thu vien Boundary Conditions
 * 1. Add Id cua dieu kien bien moi vao file bcVariables.h, vi du
 *    static constexpr int densityZeroGrad = 12;
 * 2. Them phan I/O vao DGIOLib.cpp, vi du:
 *    IO::readScalarBC::p() hoac DGIOLib.cpp/IO::readScalarBC::T() ...
 *    else if ((str0.compare("zeroGradRhoUncorrectP") == 0))
 *    { ... }
 * 3. Bo sung dieu kien bien vao ham correctPriVars() hoac correctPriVarsGrad() trong file BCSupportFncs.cpp
 * 4. Bo sung cau lenh switch sang new Id trong ham correctPressure(), correctTemperature() hoac correctVelocity()
 * VA cac ham correct grad tuong ung, correctPressureGrad(), correctTemperatureGrad() hoac correctVelocityGrad(), vi du:
 *      case BCVars::pressureBCId::zeroGradRhoUncorrectP:
        {
            //zeroGradRho
            pM = zeroGradient_scalar(pP); //==> update lai p theo TM va rhoP o ngoai ham nay
        }
        break;
*/

#ifndef DGBCSLIB_H_INCLUDED
#define DGBCSLIB_H_INCLUDED
#include <vector>
/*NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
*/

/*Function applies advective and diffusive boundary conditions to edges in NSF equation, it returns values of fluxes at Gauss point*/
std::vector<double> NSFEqBCsImplement(int element, int edge, int nG);
/*Cac ham con cua NSFEqBCsImplement*/
void NSFEqBCsForNormalBoundary(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG);

/*Function applies auxilary boundary conditions to edges in auxilary equation, it returns values of auxiraly flux at Gauss point*/
std::vector<std::vector<double>> auxEqBCsImplement(int element, int edge, int nG);
/*Cac ham con cua auxEqBCsImplement*/
void auxEqBCsForNormalBoundary(std::vector<std::vector<double>> &Fluxes, int element, int edge, int edgeGrp, int nG);

//Implement bondary condition of Rho (use when massDiffusion is on)
std::tuple<double, double> rhoBCsImplement(int edge, int nG);

#endif // DGBCSLIB_H_INCLUDED
