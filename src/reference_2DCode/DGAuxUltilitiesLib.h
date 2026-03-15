#ifndef DGAUXULTILITIESLIB_H_INCLUDED
#define DGAUXULTILITIESLIB_H_INCLUDED
#include <tuple>
#include <vector>
#include <string>

int calcArrId(int id1, int id2, int length);

namespace auxUlti
{
	/*Funtion finds order of edge respective to element*/
	int findEdgeOrder(int element, int edge);

	/*Function checks element's type, return 3 if element is tri, 4 if element is quad*/
    int checkType(int element);

	/*Function gets coordinates of index'th vertex of element*/
	std::tuple<double, double> getElemCornerCoord(int elem, int index);

	/*Function calculates primary variables from conservative variables*/
	//void ConserToPri();

	/*Function gets working directory*/
	std::string workingdir();

	/*Function return true if considering element is master of considering edge, otherwise it return false*/
	bool checkMaster(int elem, int edge);

    /**
     * @brief Function check boundary type of edge (wall, patch, symmetry, matched).
     *
     * @param edgeId: global edge Id.
     * @return Id of BC type (compare to meshVar::BCTypeID Ids to determine BC Type).
     */
    int checkBCTypeOfEdge(int edgeId);

	/*Function gets normal vector component (nx or ny) from normalVector array*/
	double getNormVectorComp(int elem, int edge, int dir);

	/*Function gets (a,b) coordinates of Gauss point on surface*/
	std::tuple<double, double> getGaussSurfCoor(int edge, int elem, int nG);

	/*Function gets (a,b) coordinates of Gauss point on surface of cell master*/
	std::tuple<double, double> getGaussSurfCoorMaster(int edge, int elem, int nG);

	/*Function executes file .exe*/
	void openFileEXE(std::string location);

	/*Function gets value of all order of accuracy of conservative variables at inputted element and returns output as a vector
	Type 1: rho
	2: rhou
	3: rhov
	4: rhoE*/
	std::vector<double> getElementConserValuesOfOrder(int element, int type);

	/*Function gets value of all order of accuracy of auxilary variables at inputted element and returns output as a vector
	Type 1: drho
	2: drhou
	3: drhov
	4: drhoE
	
	dir 1: Ox
	dir 2: Oy*/
	std::vector<double> getElementAuxValuesOfOrder(int element, int type, int dir);

	/*Function gets (a,b) coordinates of Gauss point at inside element*/
	std::tuple<double, double> getGaussCoor(int na, int nb);

	/*Function gets group index of edge*/
	int getGrpOfEdge(int edge);

	/*Function gets boundary type of edge*/
	int getBCType(int edge);

	/*Function returns true if problem is subsonic*/
	bool checkSubSonic();

	/*Function checks subsonic flow locally*/
	bool checkSubSonicLocally(double TVal, double uVal, double vVal);

    bool checkNonEqmBCAvailable();

    bool checkTimeVaryingBCAvailable();

	/*Function returns master element and servant element of edge*/
	std::tuple<int, int> getMasterServantOfEdge(int edge);

	//Function returns cell centroid coordinates and size (cell area)
	std::tuple<double, double, double> getCellMetrics(int element);

    //void resize1DArray(double*Array, int row, double initialValue);

    //void resize1DIntArray(int*Array, int row, int initialValue);

    void initialize1DArray(double*Array, int row, double initialValue);

    void initialize1DIntArray(int*Array, int row, int initialValue);

    void initialize1DBoolArray(bool*Array, int row, bool initialValue);

	/*Function resize 2D array type double*/
    double** resize2DArray(int row, int column, double initialValue);

	/*Function resize 3D array*/
    void resize3DArray(std::vector<std::vector<std::vector<double>>> &Array, int direct1, int direct2, int direct3);

	/*Function resizes all dynamic arrays, it helps to reduce amount of consumed RAM*/
	void resizeDGArrays();

	/*Function resize 2D array type int*/
     int** resize2DIntArray(int row, int column, int initialValue);

	/*Function computes coordinates of Gauss point on all edges*/
	void mappingEdges();

	//This function supports for inverse coodinates mapping
	std::vector<std::vector<double>> getVectorGaussSurfCoor(int edge, int elem);

	//Function returns location of input edge on BC values array
    int getAdressOfBCEdgesOnBCValsArray(int edge);

    //Function gets globle edge id from local BC edge id
    int getGlobalEdgeIdFromLocalBCEdgeId(int localBCEdgeId);

	//Function gets centroid coordinates of inputted cell
	std::tuple<double, double> getCellCentroid(int element);

	//Function deletes 1D vector and frees its memory
	void clear1DIntVector(std::vector<int>&vector);

    //Function deletes 2D vector and frees its memory
    void clear2DIntVector(std::vector<std::vector<int>>&vector);

	//Function returns neighbor element of input element which shares input edge with
	int getNeighborElement(int element, int edge);

	//Function gets edge of input element which has input order
	int getEdgeHasInputOrderOfElement(int element, int inputEdgeOrder);

	//Function gets value of all order of accuracy of residuals at inputted element and returns output as a vector
	std::vector<double> getResidualValuesOfOrder(int element, int type);

	//Function finds order of input point respect to input element
	int findVertexOrder(int point, int element);

	//Function add row with numCol to input array
    void addRowTo2DIntArray(std::vector<std::vector<int>> &Array, int numCol);
    void addRowTo2DDoubleArray(std::vector<std::vector<double>> &Array, int numCol);

    //Function create folder at input location
    void createFolder(std::string location, bool passExit);

    std::string createTimeStepFolder(int iter, std::string option);

    //Ham lay gia tri UMean (U order 0) cua 1 cell
    void getVectorUMeanOfCell(int element, std::vector<double> &U);

    std::tuple<double, double> getUAtInterfaces(int edge, int element, int nG, int valType);
    double getUPlusAtBC(int edge, int nG, int valType);

    std::tuple<double, double> getTAtInterfaces(int edge, int element, int nG);
    double getTPlusAtBC(int edge, int nG);

    void saveUAtBCToSurfaceFields(int edge, int nG, std::vector<double>&UPlus, std::vector<double>&UMinus);

    std::vector<double> getElementAuxValuesOfOrder_BR2_vol(int element, int type, int dir);

    std::vector<double> getElementAuxValuesOfOrder_BR2_sur(int edge, int element, int type, int dir);

    void copyFolder(std::string source, std::string destination);

    void copyFile(std::string source, std::string destination);

	//Auxilary functions support for postProcessing
	namespace postProcess
	{
		//Function gets vector of element which are surrounding inputted point
		std::vector<int> getElementsSurroundingPoint(int point);

		//Function finds index of inputted integer in inputted array
		std::tuple<double, double> findPointCoorInStandardSpace(int point, int element);
	}

    void checkInforBeforeRunning();

    void getCommand();

    void shrink2DIntVector(std::vector<std::vector<int>>&vector, int numRow);

    void releaseMemory();

    void resizeTemporaryArrays();

    std::tuple<int, bool> lookForDataOfKeyword(std::string fileLoc, std::string inputKeyWord);

    void resizeRequiredArrays();
}
#endif // DGAUXULTILITIESLIB_H_INCLUDED
