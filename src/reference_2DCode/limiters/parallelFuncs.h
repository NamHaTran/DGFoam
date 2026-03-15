#ifndef PARALLELFUNCS_H
#define PARALLELFUNCS_H
#include <tuple>
namespace limiter_parallelFuncs
{
    //void sendRecvSurfGaussPtArray(double**Var,double**Buffer);

    //void sendRecvTroubleCellGaussValue(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG);

    namespace massDiff {
        void sendRecvGaussPtValues(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG);

        void sendRecvSurfGaussArray(double**Var,double**Buffer);

        void synchCellStatus();

        void sendRecvCellStatusBtwProcs(int sendingProc, int receivingProc, bool*Var,bool*Buffer);

        void countTrbEdge();
    }
}

#endif // PARALLELFUNCS_H
