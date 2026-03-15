#include "symmetry.h"
#include <vector>
#include "DGAuxUltilitiesLib.h"
#include "BCSupportFncs.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGMath.h"
#include "bcVariables.h"

//Debug
#include "debuggingFuncs.h"
#include <iostream>

//Durst Model
#include "./extNSFEqns/FranzDurst/DurstModel.h"

// symmetry
/*
void symmetry_vector(std::vector<double> &vectorM, const std::vector<double> &vectorP,  const std::vector<double> &n)
{
    vectorM[0] = vectorP[0] - 2 * (vectorP[0] * n[0] + vectorP[1] * n[1])*n[0];
    vectorM[1] = vectorP[1] - 2 * (vectorP[0] * n[0] + vectorP[1] * n[1])*n[1];


    double normP(vectorP[0]*n[0]+vectorP[1]*n[1]),
            tangP(vectorP[0]*n[1]-vectorP[1]*n[0]);
    //Nguyen tac: normM = -normP, tangM = tangP
    vectorM[0] = -normP*n[0]+tangP*n[1];
    vectorM[1] = (-normP - vectorM[0]*n[0]);

    //Correct vectorM[1]
    //Vi vectorM = (-normP - vectorM[0]*n[0])/n[1] nen phai correct de tranh chia cho 0
    if ((fabs(n[1])<1e-5) && (fabs(vectorM[1])<1e-5))
    {
        vectorM[1] = 0;
    }
    else
    {
        vectorM[1] /=n[1];
    }
}

double symmetry_scalar(double scalarP)
{
    return (scalarP);
}*/


void symmetryBC_auxEqs(std::vector <std::vector<double>> &Fluxes, int element, int edge, int edgeGrp, int nG)
{
    std::vector<double>
        UPlus(4, 0.0),
        UMinus(4, 0.0),
        norm(2, 0.0);
    double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2)), muP;
    std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
    norm[0] = nx;
    norm[1] = ny;

    //Compute plus values--------------------------------------------------------
    BCSupportFncs::auxilaryBCs::getUPlus(edge,nG,UPlus);
    //BCSupportFncs::auxilaryBCs::calcUMean(element,UPlus);

    //Compute minus values--------------------------------------------------------
    UMinus[0] = UPlus[0];
    if (BCVars::DirichletAppMethGeneralBCStrong[edgeGrp-1])
    {
        UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
        UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
    }
    else
    {
        /*
        double tang(-UPlus[1]*ny+UPlus[2]*nx);
        UMinus[1]=-tang*ny;
        UMinus[2]=tang*nx;*/

        UMinus[1] = UPlus[1] - (UPlus[1] * nx + UPlus[2] * ny)*nx;
        UMinus[2] = UPlus[2] - (UPlus[1] * nx + UPlus[2] * ny)*ny;
    }

    UMinus[3] = UPlus[3];
    //----------------------------------------------------------------------------

    //Nhan mu vao vector flux truoc khi return
    muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);

    //Luu TMinus (TMinus = TPlus)
    surfaceFields::T[edge][nG+mathVar::nGauss1D+1]=surfaceFields::T[edge][nG];

    for (int i = 0; i < 4; i++)
    {
        Fluxes[i][0] = UPlus[i]*muP;
        Fluxes[i][1] = UMinus[i]*muP;
    }

    //Save U+ and U- at boundary to arrays
    auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
}

void symmetryBC_NSFEqs(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG)
{
    std::vector<double>
        UPlus(4, 0.0),
        UMinus(4, 0.0),
        dUXPlus(4, 0.0),
        dUYPlus(4, 0.0),
        dUXMinus(4, 0.0),
        dUYMinus(4, 0.0),
        norm(2, 0.0);
    double a(0.0), b(0.0), TPlus(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
    std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
    norm[0] = nx;
    norm[1] = ny;

    //Compute U+
    for (int i = 0; i < 4; i++)
    {
        std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
    }
    TPlus=surfaceFields::T[edge][nG];

    //Compute dU+/-
    if (flowProperties::viscous)
    {
        BCSupportFncs::NSFEqBCs::getdUPlus(edge,nG,dUXPlus,dUYPlus);
        //BCSupportFncs::NSFEqBCs::calcdUMean(edge,element,dUXPlus,dUYPlus);
    }

    //double tangGrad(0.0);

    if (BCVars::NewmannAppMethGradGeneralBCStrong[edgeGrp-1])
    {
        for (int i = 0; i < 4; i++)
        {
            dUXMinus[i] = dUXPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
            dUYMinus[i] = dUYPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;
        }
    }
    else
    {
        for (int i = 0; i < 4; i++)
        {
            /*
            tangGrad = -dUXPlus[i]*ny+dUYPlus[i]*nx;
            dUXMinus[i]=-tangGrad*ny;
            dUYMinus[i]=tangGrad*nx;*/
            dUXMinus[i] = dUXPlus[i] - (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
            dUYMinus[i] = dUYPlus[i] - (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;
        }
    }

    //DURST MODEL-------------------------------------------------------
    //Set flag isSymmetry to True to correct self-diffusive flux at symmetry BC
    extNSF_Durst::isSymmetry=true;
    //DURST MODEL-------------------------------------------------------

    BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, TPlus, TPlus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);

    //DURST MODEL-------------------------------------------------------
    //Reset flag isSymmetry
    extNSF_Durst::isSymmetry=false;
    //DURST MODEL-------------------------------------------------------
}
