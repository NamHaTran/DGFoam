#include "./limiters/parallelFuncs.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include <mpi.h>
#include "./parallelFunctions/GaussPointData.h"
#include "./parallelFunctions/generalParallelFuncs.h"
#include "DGAuxUltilitiesLib.h"

#include "./limiters/massDiffusion/massDiffusion.h"

//Debug
#include <iostream>
#include "DGIOLib.h"

namespace limiter_parallelFuncs
{
    namespace massDiff {

        /**
         * @brief Function synchs cells (which has matched BC edge) status: be trouble or not trouble.
         *
         * Function also check whether the process has any trouble cells at "matched" BC. If no, parallel functions will not be ran.
         *
         * @param sendingProc: id of sending proc.
         * @param receivingProc: id of receiving proc.
         * @param Var: markerOfTrbCellAtMatchedBC array.
         * @param Buffer: markerOfTrbCellAtMatchedBC array.
         */
        void sendRecvCellStatusBtwProcs(int sendingProc, int receivingProc, bool*Var,bool*Buffer)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            //Declare dynamic array
            MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

            int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0), edgeId;

            if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
            {
                //Send
                for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
                {
                    destination=meshVar::meshConnection[iBCedge][1];
                    neighborCellId=meshVar::meshConnection[iBCedge][2];

                    if (destination==receivingProc)
                    {
                        //define tags
                        tag_sent = neighborCellId;

                        MPI_Isend(&Var[iBCedge], 1, MPI_CXX_BOOL, destination, tag_sent, MPI_COMM_WORLD, &request[requestCount]);
                        requestCount++;
                    }
                }
                //Wait for all request fullfilled
                MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
            }
            else if (systemVar::currentProc==receivingProc)
            {
                bool recvVal(false);

                //Receive
                for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++) {
                    source=meshVar::meshConnection[iBCedge][1];
                    cellId=meshVar::meshConnection[iBCedge][0];

                    if (source==sendingProc)
                    {
                        //define tags
                        tag_recv = cellId;

                        MPI_Irecv(&Buffer[iBCedge], 1, MPI_CXX_BOOL, source, tag_recv, MPI_COMM_WORLD, &request[requestCount]);
                        requestCount++;
                    }
                }
                //Wait for all request fullfilled
                MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            delete [] request;
        }

        void sendRecvGaussPtValues(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            //Declare dynamic array

            if (limiter::massDiffusion::trbCellAtMatchedBC)
            {
                MPI_Request *request = new MPI_Request[limiter::massDiffusion::numOfTrbEdge];

                int requestCount(0);

                if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
                {
                    //Send
                    for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
                    {
                        parallelFuncs_GaussPt::sendDataOneEdge(iBCedge,receivingProc,Var,nG,request,requestCount);
                    }
                    //Wait for all request fullfilled
                    MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
                }
                else if (systemVar::currentProc==receivingProc)
                {
                    //Receive
                    for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++) {
                        parallelFuncs_GaussPt::recvDataOneEdge(iBCedge,sendingProc,Buffer,mathVar::nGauss1D + nG + 1,request,requestCount);
                    }
                    //Wait for all request fullfilled
                    MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
                }
                delete [] request;
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

        void sendRecvSurfGaussArray(double**Var,double**Buffer)
        {
            int sendingProc, receivingProc;
            for (int i=0; i<systemVar::sendRecvOrder_length; i++)
            {
                sendingProc=systemVar::sendRecvOrder[i][0];
                receivingProc=systemVar::sendRecvOrder[i][1];
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    limiter_parallelFuncs::massDiff::sendRecvGaussPtValues(sendingProc,receivingProc,Var,Buffer,nG);
                }
            }
        }

        /**
         * @brief Function synchs status between all processors.
         */
        void synchCellStatus()
        {
            int sendingProc, receivingProc;
            for (int i=0; i<systemVar::sendRecvOrder_length; i++)
            {
                sendingProc=systemVar::sendRecvOrder[i][0];
                receivingProc=systemVar::sendRecvOrder[i][1];

                limiter_parallelFuncs::massDiff::sendRecvCellStatusBtwProcs(
                            sendingProc,receivingProc,
                            limiter::massDiffusion::markerOfTrbCellAtMatchedBC,
                            limiter::massDiffusion::markerOfTrbCellAtMatchedBC_buffer);
            }

            //Synch
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
            {
                limiter::massDiffusion::markerOfTrbCellAtMatchedBC[iBCedge]=
                        limiter::massDiffusion::markerOfTrbCellAtMatchedBC[iBCedge]
                        ||
                        limiter::massDiffusion::markerOfTrbCellAtMatchedBC_buffer[iBCedge];
            }

            //Count trouble edges after synched
            limiter_parallelFuncs::massDiff::countTrbEdge();
        }

        /**
         * @brief Function counts number of trouble edges on "matched" BC after synching status between all processors.
         */
        void countTrbEdge()
        {
            //Tao bien temporary numTrbEdge chu khong count truc tiep vao limiter::massDiffusion::numOfTrbEdge de reset lai num of trb edge sau moi lan chay ham countTrbEdge.
            int numTrbEdge(0);

            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
            {
                if (limiter::massDiffusion::markerOfTrbCellAtMatchedBC[iBCedge])
                    numTrbEdge++;
            }

            limiter::massDiffusion::numOfTrbEdge=numTrbEdge;
            if (numTrbEdge>0)
                limiter::massDiffusion::trbCellAtMatchedBC=true;
            else
                limiter::massDiffusion::trbCellAtMatchedBC=false;

            //IO::write1DBoolVectorToFile("/mnt/d/Kei_Nam/Documents/DGSolver/DG2DSolver_newSendRecv/output", "file"+std::to_string(systemVar::currentProc),limiter::massDiffusion::markerOfTrbCellAtMatchedBC,meshVar::numBCEdges); //DEBUG_COMM
        }
    }
}
