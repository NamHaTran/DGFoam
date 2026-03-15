#ifndef DGMESHREADERLIB_H_INCLUDED
#define DGMESHREADERLIB_H_INCLUDED
#include <tuple>  //Include this for returning multiple values in function
#include <vector>
#include <string>
namespace MshReader
{
	/*Function creates matrixes which save mesh data get from files, by using EleSurPtProcess()*/
	void meshProcess();

	/*Function gets informations of elements surrounding points*/
	void EleSurPt();

	/*Function gets informations of points surrounding points*/
	void PtsSurPt();

	/*Function gets informations of elements surrounding element*/
	void ElemsSurElem();

	/*Function gets informations of edges*/
	void EdgesInfor();
	void EdgesOfElem();

	/*Calculate normal vector of each face (edge)*/
	void GetNormalVector();

	/*Function supports for EdgesInfor(), it returns index of checking number in input iarray*/
	int findIndex(int number, int iarray[], int size);

	//Note: run this function AFTER mesh processing is DONE
	void getBoundaryPoints();

	/*Child functions*/
	int CheckConnection(int point, int helpArray[], int length);
	int checkIndividualEdge(int rootPt, int tipPt);
	int findIndex(int number, int iarray[], int size);
	void getBcGrpTp(int ipoin, int jpoin, int ninpoed);
	std::tuple<double, double> calcNormVector(int point1, int point1Indice, int point2, int point2Indice, int type);

	void sortPointsOfElements();

    void getMidpointOfBCEdge(std::string mode);

    void fixCellVerticesOrder();
}

namespace MshExporter
{
    void exportMeshToMetis();

    void testMeshPartitionResult();
}

namespace decomposePart {
    std::vector<int> loadPartitionedMesh();

    std::vector<int> findLocalIdOfPts();

    void findEdgeWithBCTypeMatched();

    std::vector<int> getMeshInforOfRanks(std::vector<std::vector<std::vector<double>>>&Points, std::vector<std::vector<std::vector<int>>>&Elem1D, std::vector<std::vector<std::vector<int>>>&Elem2D, std::vector<std::vector<std::vector<int>>>&meshConnection);

    void decomposingMesh();

    void exportPartitionedMesh(int rank, int npoin, int nelem2D, std::vector<std::vector<double>>&Points, std::vector<std::vector<int>>&Elements2D);

    void decomposingTime0(std::string Loc);

    void decomposingLatestTime();

    void distributingDiscretedVar();

    void fixPartitionedMesh();

    //void readLocalIdAndRankOf2DElem(int *rankOf2DElem,int *Elem2DlocalIdWithRank);
}

#endif // DGMESHREADERLIB_H_INCLUDED
