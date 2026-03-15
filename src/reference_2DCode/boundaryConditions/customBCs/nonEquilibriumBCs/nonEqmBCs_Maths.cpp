#include "nonEqmBCs_Maths.h"
#include "nonEqmBCs_Vars.h"
#include "DGMath.h"
#include "DGAuxUltilitiesLib.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"

#include <iostream>

namespace geometricOp {
    void findGaussPtsOnWallParameters()
    {
        int globleEdge(0), element(0);
        double a, b, xC, yC, xGauss, yGauss;
        for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
        {
            globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
            int bcType(auxUlti::checkBCTypeOfEdge(globleEdge));
            std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(globleEdge);
            if (bcType==meshVar::BCTypeID::wall) //type wall
            {
                for (int nG=0; nG<=mathVar::nGauss; nG++)
                {
                    //Get Gauss point's coordinates
                    std::tie(a, b) = auxUlti::getGaussSurfCoor(globleEdge, element, nG);

                    //Map from standard to real coor sys
                    std::tie(xGauss, yGauss) = math::directMapping(element, a, b);

                    meshVar::GaussPtsOnBCEdge_x[ilocalEdge][nG] = xGauss;
                    meshVar::GaussPtsOnBCEdge_y[ilocalEdge][nG] = yGauss;

                    //Get centroid
                    xC=meshVar::geoCenter[element][0];
                    yC=meshVar::geoCenter[element][1];

                    //Calc distance from centroid to Gauss point
                    meshVar::distanceFromGaussPtsToCentroid[ilocalEdge][nG] = math::geometricOp::calDistBetween2Points(xC,yC,xGauss,yGauss);

                    //Find unit vector (Centroid)---------->(Gauss on edge)
                    std::tie(meshVar::GaussPtsOnBCEdge_unitVector_x[ilocalEdge][nG],
                             meshVar::GaussPtsOnBCEdge_unitVector_y[ilocalEdge][nG])
                            =math::geometricOp::calcUnitVectorOf2Points(xC,yC,xGauss,yGauss);
                }
            }
        }
    }
}
