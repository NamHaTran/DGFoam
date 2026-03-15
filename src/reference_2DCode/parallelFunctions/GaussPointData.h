#ifndef GAUSSPOINTDATA_H
#define GAUSSPOINTDATA_H
#include <vector>
#include <string>
#include <mpi.h>

namespace parallelFuncs_GaussPt
{
    void sendRecvMatchedBCGaussPtValues(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG);

    void sendRecvMatchedBCGaussPtCoors(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG);

    void sendRecvSurfGaussArray(double**Var,double**Buffer);

    void synchSurfGaussU();

    void synchGaussPtCoors();

    void sendDataOneEdge(int BCEdgeId, int receivingProc, double**Var, int adress2, MPI_Request*request, int &requestCount);

    void recvDataOneEdge(int BCEdgeId, int sendingProc, double**Var, int adress2, MPI_Request*request, int &requestCount);
}
#endif // GAUSSPOINTDATA_H
