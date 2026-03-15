#include "generalParallelFuncs.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <mpi.h>
#include "parallelVariables.h"

namespace parallelFuncs_Gen
{
    /**
     * @brief Function get Gauss point coordinates of neighbor cell which belongs to another master core.
     * @param loc: adress of neighor cell in buffer arrays, can be got using parallelFunctions::getGaussPointCoorsOfNeighborCell function.
     * @param nG: Gauss point Id.
     * @return
     */
    /*
    std::tuple<double,double> getGaussPointCoorsOfNeighborCell(int loc, int nG)
    {
        double a(parallelBuffer::aCoor[loc][nG]), b(parallelBuffer::bCoor[loc][nG]);
        return std::make_tuple(a,b);
    }*/

    /**
     * @brief Function prepare case for parallel running
     */
    void prepareParallelCase()
    {
        MPI_Init(NULL, NULL);
        int maxNode;
        MPI_Comm_size(MPI_COMM_WORLD, &maxNode);
        if (maxNode<systemVar::totalProc)
        {
            std::cout<<"Number of processors ("<<systemVar::totalProc<<") exceeds maximum number of available processors of system ("<<maxNode<<").\n";
            std::cout << "DGSolver will exit after you hit return.\n";
            exit(EXIT_FAILURE);
        }
        MPI_Comm_rank(MPI_COMM_WORLD, &systemVar::currentProc);
    }

    void sendReceiveMeshData(int vertex, int dir, double**Buffer)
    {
        //dir = 1: xCoor
        //dir = 2: yCoor

        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];

            MPI_Barrier(MPI_COMM_WORLD);
            //Declare dynamic array
            MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

            int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0);
            if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
            {
                //Send
                for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
                {
                    destination=meshVar::meshConnection[iBCedge][1];
                    cellId=meshVar::meshConnection[iBCedge][0];
                    neighborCellId=meshVar::meshConnection[iBCedge][2];

                    int elemType(auxUlti::checkType(cellId)), ptId(0);

                    if (destination==receivingProc)
                    {
                        tag_sent=neighborCellId;
                        //send toa do dinh cua cell
                        if (vertex<elemType)
                        {
                            ptId=meshVar::Elements2D[cellId][vertex];
                            MPI_Isend(&meshVar::Points[ptId][dir-1], 1, MPI_DOUBLE, destination, tag_sent*100+10+vertex, MPI_COMM_WORLD, &request[requestCount]);
                            requestCount++;
                        }
                    }
                }
                //Wait for all request fullfilled
                MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
            }
            else if (systemVar::currentProc==receivingProc)
            {
                //Receive
                for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++) {
                    source=meshVar::meshConnection[iBCedge][1];
                    cellId=meshVar::meshConnection[iBCedge][0];
                    //neighborCellId=meshVar::meshConnection[iBCedge][2];

                    //tag = neighborCellId*10 + destination
                    if (source==sendingProc)
                    {
                        int elemType(auxUlti::checkType(cellId));

                        tag_recv=cellId;
                        //receive toa do dinh cua cell
                        if (vertex<elemType)
                        {
                            MPI_Irecv(&Buffer[iBCedge][vertex], 1, MPI_DOUBLE, source, tag_recv*100+10+vertex, MPI_COMM_WORLD, &request[requestCount]);
                            requestCount++;
                        }
                    }
                }
                //Wait for all request fullfilled
                MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            delete [] request;
        }
    }

    void sendString(std::string content, int destination, int tag)
    {
        MPI_Send(content.c_str(), content.size(), MPI_CHAR, destination, tag, MPI_COMM_WORLD);
    }

    std::string receiveString(int source, int tag)
    {
        //From internet (giu nguyen khong sua)
        /*
        MPI_Status status;
        MPI_Probe(source, tag, MPI_COMM_WORLD, &status);
        int l;
        MPI_Get_count(&status, MPI_INT, &l);
        char *buf = new char[l];
        MPI_Recv(buf, l, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
        std::string str(buf, l);
        delete[] buf;
        return str;*/
        //From internet (giu nguyen khong sua)
        MPI::Status status;
        MPI::COMM_WORLD.Probe(source, tag, status);
        int l = status.Get_count(MPI_CHAR);
        char *buf = new char[l];
        MPI::COMM_WORLD.Recv(buf, l, MPI_CHAR, source, tag, status);
        std::string str(buf, l);
        delete[] buf;
        return str;
    }

    void resizeMeshParallelBuffers()
    {
        //Surface Gauss point coordinates
        //parallelBuffer::aCoor= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
        //parallelBuffer::bCoor= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
        //Neighbor cell vertexes
        parallelBuffer::xCoor= auxUlti::resize2DArray(meshVar::numBCEdges, 4,0.0);
        parallelBuffer::yCoor= auxUlti::resize2DArray(meshVar::numBCEdges, 4,0.0);

        parallelBuffer::elemType = new int [meshVar::numBCEdges];
        auxUlti::initialize1DIntArray(parallelBuffer::elemType,meshVar::numBCEdges,3);
    }

    /*
    void resizeGaussPtParallelBuffers()
    {
        //Buffer (parallel computing)
        //Conservative variables
        parallelBuffer::surfaceGaussPt::rho= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,iniValues::rhoIni);
        parallelBuffer::surfaceGaussPt::rhou= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
        parallelBuffer::surfaceGaussPt::rhov= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
        parallelBuffer::surfaceGaussPt::rhoE= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);

        //Temperature
        parallelBuffer::surfaceGaussPt::T= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,iniValues::TIni);

        //Auxilary variables
        if (flowProperties::viscous)
        {
            parallelBuffer::surfaceGaussPt::drhoX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
            parallelBuffer::surfaceGaussPt::drhouX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
            parallelBuffer::surfaceGaussPt::drhovX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
            parallelBuffer::surfaceGaussPt::drhoEX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
            parallelBuffer::surfaceGaussPt::drhoY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
            parallelBuffer::surfaceGaussPt::drhouY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
            parallelBuffer::surfaceGaussPt::drhovY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
            parallelBuffer::surfaceGaussPt::drhoEY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
        }
    }*/

    int getMatchedEdgeIdOfCell(int element)
    {
        int elemType(auxUlti::checkType(element)), edgeId(-1);
        for (int edgeOrder=0; edgeOrder<elemType; edgeOrder++)
        {
            edgeId = meshVar::inedel[element][edgeOrder];
            if (auxUlti::getBCType(edgeId)==4)
            {
                break;
            }
        }
        return edgeId;
    }

    /**
     * @brief Function calculates sum of int variable over all processors and send back result to all processors.
     * @param var: variable needed to sum.
     * @return
     */
    int sumIntOverProcs(int var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        int result(0);
        if (systemVar::currentProc==0)
        {
            result=var;
            int recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_INT, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                result+=recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_INT, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return result;
    }

    /**
     * @brief Function finds min of int variable over all processors and send back result to all processors.
     * @param var: variable needed to find min.
     * @return
     */
    int minIntOverProcs(int var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        int result(1e9);
        if (systemVar::currentProc==0)
        {
            result=var;
            int recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_INT, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Operator
                if (recv<result)
                    result=recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_INT, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return result;
    }

    /**
     * @brief Function finds max of int variable over all processors and send back result to all processors.
     * @param var: variable needed to find max.
     * @return
     */
    int maxIntOverProcs(int var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        int result(-1e9);
        if (systemVar::currentProc==0)
        {
            result=var;
            int recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_INT, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Operator
                if (recv>result)
                    result=recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_INT, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return result;
    }

    /**
     * @brief Function calculates sum of double variable over all processors and send back result to all processors.
     * @param var: variable needed to sum.
     * @return
     */
    double sumDoubleOverProcs(double var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        double result(0);
        if (systemVar::currentProc==0)
        {
            result=var;
            double recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_DOUBLE, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                result=recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_DOUBLE, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return result;
    }

    /**
     * @brief Function finds min of double variable over all processors and send back result to all processors.
     * @param var: variable needed to find min.
     * @return
     */
    double minDoubleOverProcs(double var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        double result(1e9);
        if (systemVar::currentProc==0)
        {
            result=var;
            double recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_DOUBLE, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Operator
                if (recv<result)
                    result=recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_DOUBLE, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return result;
    }

    /**
     * @brief Function finds max of double variable over all processors and send back result to all processors.
     * @param var: variable needed to find max.
     * @return
     */
    double maxDoubleOverProcs(double var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        double result(1e9);
        if (systemVar::currentProc==0)
        {
            result=var;
            double recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_DOUBLE, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Operator
                if (recv>result)
                    result=recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_DOUBLE, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return result;
    }

    /**
     * @brief Function calculates 'OR' operator of variable over all processors and send back result to all processors.
     * @param var: variable needed to do 'OR'.
     * @return
     */
    bool bool_OR_OverProcs(bool var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        bool result(false);
        if (systemVar::currentProc==0)
        {
            result=var;
            bool recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_CXX_BOOL, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Operator
                result=result||recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_CXX_BOOL, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        return result;
    }

    /**
     * @brief Function calculates 'AND' operator of variable over all processors and send back result to all processors.
     * @param var: variable needed to do 'AND'.
     * @return
     */
    bool bool_AND_OverProcs(bool var)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        bool result(false);
        if (systemVar::currentProc==0)
        {
            result=var;
            bool recv;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Recv(&recv, 1, MPI_CXX_BOOL, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Operator
                result=result&&recv;
            }
        }
        else {
            MPI_Send(&var, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD);
        }

        //send back result to all processors
        if (systemVar::currentProc==0)
        {
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                MPI_Send(&result, 1, MPI_CXX_BOOL, irank, 0, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(&result, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        return result;
    }
}
