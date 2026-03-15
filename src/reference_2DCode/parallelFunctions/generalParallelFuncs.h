#ifndef GENERALPARALLELFUNCS_H
#define GENERALPARALLELFUNCS_H
#include <vector>
#include <string>
namespace parallelFuncs_Gen
{
    //std::tuple<double,double> getGaussPointCoorsOfNeighborCell(int loc, int nG);
    void prepareParallelCase();
    void sendString(std::string content, int destination, int tag);
    std::string receiveString(int source, int tag);
    void sendReceiveMeshData(int vertex, int dir, double**Buffer);
    int getMatchedEdgeIdOfCell(int element);

    void resizeMeshParallelBuffers();

    //void resizeGaussPtParallelBuffers();

    int sumIntOverProcs(int var);

    int maxIntOverProcs(int var);

    int minIntOverProcs(int var);

    double sumDoubleOverProcs(double var);

    double maxDoubleOverProcs(double var);

    double minDoubleOverProcs(double var);

    bool bool_OR_OverProcs(bool var);

    bool bool_AND_OverProcs(bool var);
}
#endif // GENERALPARALLELFUNCS_H
