#include "matched.h"
#include <vector>
#include "DGAuxUltilitiesLib.h"
#include "BCSupportFncs.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGMath.h"

void matchedBC_auxEqs(std::vector <std::vector<double>> &Fluxes, int edge, int nG)
{
    std::vector<double>
        UPlus(4, 0.0),
        UMinus(4, 0.0);
    double muP, muM;

    //Compute plus values--------------------------------------------------------
    BCSupportFncs::auxilaryBCs::getUPlus(edge,nG,UPlus);

    BCSupportFncs::parallel::getUMinus_convective(edge,nG,UMinus);

    //Nhan mu vao vector flux truoc khi return
    muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
    muM=math::CalcVisCoef(surfaceFields::T[edge][nG+mathVar::nGauss1D+1]);
    for (int i = 0; i < 4; i++)
    {
        Fluxes[i][0] = UPlus[i]*muP;
        Fluxes[i][1] = UMinus[i]*muM;
    }
}

void matchedBC_NSFEqs(std::vector<double> &Fluxes, int element, int edge, int nG)
{
    std::vector<double>
        UPlus(4, 0.0),
        UMinus(4, 0.0),
        dUXPlus(4, 0.0),
        dUYPlus(4, 0.0),
        dUXMinus(4, 0.0),
        dUYMinus(4, 0.0),
        norm(2, 0.0);
    double TPlus(surfaceFields::T[edge][nG]), TMinus(surfaceFields::T[edge][nG+mathVar::nGauss1D+1]), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
    norm[0] = nx;
    norm[1] = ny;

    //Compute U+
    for (int i = 0; i < 4; i++)
    {
        std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
    }

    //Compute dU+/-
    if (flowProperties::viscous)
    {
        BCSupportFncs::NSFEqBCs::getdUPlus(edge,nG,dUXPlus,dUYPlus);
        BCSupportFncs::parallel::getdUMinus(edge,nG,dUXMinus,1);
        BCSupportFncs::parallel::getdUMinus(edge,nG,dUYMinus,2);
    }

    BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
}
