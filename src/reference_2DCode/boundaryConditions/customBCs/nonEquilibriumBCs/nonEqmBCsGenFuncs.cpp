#include "nonEqmBCsGenFuncs.h"
#include "DGAuxUltilitiesLib.h"
#include "VarDeclaration.h"
#include "DGIOLib.h"
#include "dynamicVarDeclaration.h"
#include <iostream>

#include "../../bcVariables.h"

//Time-varying boundary conditions
#include "MaxwellSlip/u_MaxwellSlip.h"
#include "SmoluchowskyTJump/T_SmoluchowskyTJump.h"

namespace nonEquilibriumBCs {
    void readSurfaceValues(std::string Loc)
    {
        if (auxUlti::checkTimeVaryingBCAvailable())
        {
            std::string fileLoc(""), fileName("");

            //Read file TSurface & USurface
            fileName = "TSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            std::tie(controlFlag::fileAvailFlags::fileTSurface,SurfaceBCFields::TBc,std::ignore)=IO::read1DArray(fileLoc,fileName,false);
            fileName = "uSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            std::tie(controlFlag::fileAvailFlags::fileuSurface,SurfaceBCFields::uBc,std::ignore)=IO::read1DArray(fileLoc,fileName,false);
            fileName = "vSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            std::tie(controlFlag::fileAvailFlags::filevSurface,SurfaceBCFields::vBc,std::ignore)=IO::read1DArray(fileLoc,fileName,false);
        }
    }

    void writeSurfaceValues(std::string Loc)
    {
        if (auxUlti::checkTimeVaryingBCAvailable())
        {
            std::string fileLoc(""), fileName("");

            //Read file TSurface & USurface
            fileName = "TSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            IO::write1DDoubleVectorToFile(Loc,fileName,SurfaceBCFields::TBc,meshVar::numBCEdges);
            fileName = "uSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            IO::write1DDoubleVectorToFile(Loc,fileName,SurfaceBCFields::uBc,meshVar::numBCEdges);
            fileName = "vSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            IO::write1DDoubleVectorToFile(Loc,fileName,SurfaceBCFields::vBc,meshVar::numBCEdges);
        }
    }

    void updateBCs()
    {
        /*
         * Ham update gia tri tren cac field cua surfaceBCFields, chi dung cho cac BC bien thien theo thoi gian.
         * Hien tai, ham su dung cho temperatureJump va slip conditions
        */
        if (auxUlti::checkTimeVaryingBCAvailable() && flowProperties::viscous)
        {
            if (systemVar::currentProc==0)
            {
                std::cout<<"Updating time varying non equilibrium BCs.\n";
            }

            int globleEdge(0);
            for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
            {
                globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
                int edgeGrp(auxUlti::getGrpOfEdge(globleEdge));
                int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]);


                //Update Temperature
                if (TType == BCVars::temperatureBCId::SmoluchowskyTJump)
                {
                    //SmoluchowskyTJump::calcTJump_DGTypeExplicit(globleEdge,edgeGrp);
                    SmoluchowskyTJump::calcTJump_FDMTypeImplicit(globleEdge,edgeGrp);
                }

                //Update Velocity
                if (UType == BCVars::velocityBCId::MaxwellSlip)
                {
                    //MaxwellSlip::calcUSlip_DGTypeExplicit(globleEdge,edgeGrp);
                    MaxwellSlip::calcUSlip_FDMTypeImplicit(globleEdge,edgeGrp);
                }
            }
        }
    }
}
