#include "GaussPointData.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <mpi.h>
#include "parallelVariables.h"
#include "generalParallelFuncs.h"

#include "debuggingFuncs.h"

namespace parallelFuncs_GaussPt
{
    void sendRecvMatchedBCGaussPtCoors(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG)
    {
        /* Note: ham nay khac ham sendRecvMatchedBCGaussPtValues o cho:
         * - Send: Id cua nG la mathVar::nGauss + nG + 1 (nguoc lai ham sendRecvMatchedBCGaussPtValues)
         * - Recv: Id cua nG la nG
        */

        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

        int requestCount(0);

        if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
        {
            //Send
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
            {
                parallelFuncs_GaussPt::sendDataOneEdge(iBCedge,receivingProc,Var,mathVar::nGauss1D + nG + 1,request,requestCount);
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }
        else if (systemVar::currentProc==receivingProc)
        {
            //Receive
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++) {
                parallelFuncs_GaussPt::recvDataOneEdge(iBCedge,sendingProc,Buffer,nG,request,requestCount);
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    void sendRecvMatchedBCGaussPtValues(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

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

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    /**
     * @brief Function send/recv value at Gauss points on "matched" BC.
     */
    void sendRecvSurfGaussArray(double**Var,double**Buffer)
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];
            for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
            {
                sendRecvMatchedBCGaussPtValues(sendingProc,receivingProc,Var,Buffer,nG);
            }
        }
    }

    /**
     * @brief Function send/recv U values at Gauss points on "matched" BC.
     */
    void synchSurfGaussU()
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];

            for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
            {
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rho,surfaceFields::rho,nG);
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rhou,surfaceFields::rhou,nG);
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rhov,surfaceFields::rhov,nG);
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rhoE,surfaceFields::rhoE,nG);
            }
        }
    }

    void synchGaussPtCoors()
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];

            if (sendingProc>receivingProc)
            {
                for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
                {
                    sendRecvMatchedBCGaussPtCoors(sendingProc,receivingProc,meshVar::edgeGaussPoints_a,meshVar::edgeGaussPoints_a,nG);
                    sendRecvMatchedBCGaussPtCoors(sendingProc,receivingProc,meshVar::edgeGaussPoints_b,meshVar::edgeGaussPoints_b,nG);
                }
            }
        }
    }

    /**
     * @brief Function sends data on 1 BC edge.
     * @param BCEdgeId: Id of BC edge.
     * @param receivingProc: Id of receiving processor (which data sent from current processor is received).
     * @param Var: Array needs to send.
     * @param adress2: 2nd adress in Var array.
     * @param request: Array of request.
     * @param requestCount: Variable for counting number of receiving request (passed by reference).
     */
    void sendDataOneEdge(int BCEdgeId, int receivingProc, double**Var, int adress2, MPI_Request*request, int &requestCount)
    {
        int destination=meshVar::meshConnection[BCEdgeId][1], neighborCellId=meshVar::meshConnection[BCEdgeId][2];
        int edgeId, tag_sent;

        //tag = neighborCellId*10 + destination
        if (destination==receivingProc)
        {
            edgeId = auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(BCEdgeId);

            //define tags
            tag_sent = neighborCellId;

            MPI_Isend(&Var[edgeId][adress2], 1, MPI_DOUBLE, destination, tag_sent, MPI_COMM_WORLD, &request[requestCount]);
            requestCount++;
        }
    }

    /**
     * @brief Function receives data on 1 BC edge.
     * @param BCEdgeId: Id of BC edge.
     * @param sendingProc: Id of sending processor (from which data is received).
     * @param Var: Array which received data is saved to.
     * @param adress2: 2nd adress in Var array.
     * @param request: Array of request.
     * @param requestCount: Variable for counting number of receiving request (passed by reference).
     */
    void recvDataOneEdge(int BCEdgeId, int sendingProc, double**Var, int adress2, MPI_Request*request, int &requestCount)
    {
        int source=meshVar::meshConnection[BCEdgeId][1], cellId=meshVar::meshConnection[BCEdgeId][0];
        int edgeId, tag_recv;

        //tag = neighborCellId*10 + destination
        if (source==sendingProc)
        {
            edgeId = auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(BCEdgeId);

            //define tags
            tag_recv = cellId;

            MPI_Irecv(&Var[edgeId][adress2], 1, MPI_DOUBLE, source, tag_recv, MPI_COMM_WORLD, &request[requestCount]);
            requestCount++;
        }
    }
}
