#include "cellData.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <mpi.h>

#include "generalParallelFuncs.h"

namespace parallelFuncs_cell
{
    void sendRecvDiscretedVar(int sendingProc, int receivingProc, double**Var,double**Buffer,int order)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0);

        if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
        {
            //Send
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge)
            {
                destination=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (destination==receivingProc)
                {
                    //define tags
                    tag_sent = neighborCellId;

                    MPI_Isend(&Var[cellId][order], 1, MPI_DOUBLE, destination, tag_sent*1000+100+order, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }
        else if (systemVar::currentProc==receivingProc)
        {
            //Receive
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
                source=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (source==sendingProc)
                {
                    //define tags
                    tag_recv = cellId;

                    MPI_Irecv(&Buffer[iBCedge][order], 1, MPI_DOUBLE, source, tag_recv*1000+100+order, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    void sendRecvTheta(int sendingProc, int receivingProc, double*thetaArray, double*Buffer)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0);

        if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
        {
            //Send
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge)
            {
                destination=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (destination==receivingProc)
                {
                    //define tags
                    tag_sent = neighborCellId;

                    MPI_Isend(&thetaArray[cellId], 1, MPI_DOUBLE, destination, tag_sent*1000+140, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }
        else if (systemVar::currentProc==receivingProc)
        {
            //Receive
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
                source=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (source==sendingProc)
                {
                    //define tags
                    tag_recv = cellId;

                    MPI_Irecv(&Buffer[iBCedge], 1, MPI_DOUBLE, source, tag_recv*1000+140, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    void sendReceiveU()
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];
            for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder)
            {
                //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, rho,parallelBuffer::rho,iorder);
                //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, rhou,parallelBuffer::rhou,iorder);
                //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, rhov,parallelBuffer::rhov,iorder);
                //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, rhoE,parallelBuffer::rhoE,iorder);
            }
        }
    }

    void sendReceivedU()
    {
        if (systemVar::auxVariables==1)
        {
            int sendingProc, receivingProc;
            for (int i=0; i<systemVar::sendRecvOrder_length; i++)
            {
                sendingProc=systemVar::sendRecvOrder[i][0];
                receivingProc=systemVar::sendRecvOrder[i][1];
                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder)
                {
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoX,parallelBuffer::drhoX,iorder);
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoY,parallelBuffer::drhoY,iorder);

                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhouX,parallelBuffer::drhouX,iorder);
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhouY,parallelBuffer::drhouY,iorder);

                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhovX,parallelBuffer::drhovX,iorder);
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhovY,parallelBuffer::drhovY,iorder);

                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoEX,parallelBuffer::drhoEX,iorder);
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoEY,parallelBuffer::drhoEY,iorder);
                }
            }
        }
        else if (systemVar::auxVariables==2)
        {
            //Chua lam cho method BR2
        }
    }

    void sendReceivedRho()
    {
        if (systemVar::auxVariables==1)
        {
            int sendingProc, receivingProc;
            for (int i=0; i<systemVar::sendRecvOrder_length; i++)
            {
                sendingProc=systemVar::sendRecvOrder[i][0];
                receivingProc=systemVar::sendRecvOrder[i][1];
                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder)
                {
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoX,parallelBuffer::drhoX,iorder);
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoY,parallelBuffer::drhoY,iorder);
                }
            }
        }
        else if (systemVar::auxVariables==2)
        {
            //Chua lam cho method BR2
        }
    }

    void sendReceivedRho_forMassDiffOnly()
    {
        if (systemVar::auxVariables==1)
        {
            int sendingProc, receivingProc;
            for (int i=0; i<systemVar::sendRecvOrder_length; i++)
            {
                sendingProc=systemVar::sendRecvOrder[i][0];
                receivingProc=systemVar::sendRecvOrder[i][1];
                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder)
                {
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::massDiffusion::rhoX,parallelBuffer::massDiffusion::drhoX,iorder);
                    //parallelFuncs_cell::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::massDiffusion::rhoY,parallelBuffer::massDiffusion::drhoY,iorder);
                }
            }
        }
        else if (systemVar::auxVariables==2)
        {
            //Chua lam cho method BR2
        }
    }
}
