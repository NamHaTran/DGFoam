#include "timeVaryingBCs_GenFuncs.h"
#include "DGAuxUltilitiesLib.h"
#include "VarDeclaration.h"
#include "DGIOLib.h"
#include "dynamicVarDeclaration.h"
#include <iostream>
#include <sstream>

#include "../../bcVariables.h"
#include "./waveTransmissive/vector_waveTransmissive.h"
#include "./waveTransmissive/scalar_waveTransmissive.h"

namespace timeVaryingBCs {
    void updateBCs()
    {
        /*
         * Ham update gia tri tren cac field cua surfaceBCFields, chi dung cho cac BC bien thien theo thoi gian.
        */
        if (auxUlti::checkTimeVaryingBCAvailable())
        {
            if (systemVar::currentProc==0)
            {
                std::cout<<"Updating time varying BCs.\n";
            }

            int globleEdge(0);
            for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
            {
                globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
                int edgeGrp(auxUlti::getGrpOfEdge(globleEdge));
                int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);

                for (int nG=0; nG<=mathVar::nGauss1D; nG++)
                {
//------------------//Wave transmissive BC//----------------------------------------------------
                    //Update Velocity
                    if (UType == BCVars::velocityBCId::waveTransmissive)
                    {
                        waveTransmissive::solveVectorEqn(globleEdge,edgeGrp,nG);
                    }
                    //Update Temperature
                    if (TType == BCVars::temperatureBCId::waveTransmissive)
                    {
                        waveTransmissive::solveScalarEqn(globleEdge,edgeGrp,nG,"T");
                    }
                    //Update Pressure
                    if (pType == BCVars::pressureBCId::waveTransmissive)
                    {
                        waveTransmissive::solveScalarEqn(globleEdge,edgeGrp,nG,"p");
                    }
//----------------------------------------------------------------------------------------------
                }
            }
        }
    }
}
