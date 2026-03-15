#ifndef CELLDATA_H
#define CELLDATA_H
#include <vector>
#include <string>
#include "generalParallelFuncs.h"
namespace parallelFuncs_cell
{
    void sendRecvDiscretedVar(int sendingProc, int receivingProc, double**Var,double**Buffer,int order);
    void sendRecvTheta(int sendingProc, int receivingProc, double*thetaArray,double*Buffer);
    void sendReceiveU();
    void sendReceivedU();
    void sendReceivedRho();
    void sendReceivedRho_forMassDiffOnly();
}
#endif // CELLDATA_H
