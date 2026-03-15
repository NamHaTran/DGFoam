#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include <vector>

namespace meshVar
{
    /*Edges informations*/
    int **inpoed = new int*[1];
    /* - Cot 0 va 1: 2 point cua edge
     * - Cot 3: Kieu dieu kien bien cua edge (1: wall. 2: intFlow. 3: outFlow. 4: symmetry)
     * - Cot 4: id cua edge trong vector Element1D (Element 1D chi bao gom cac edge tai boundary)
    */

    /*Elements surrounding element*/
    int **esuel = new int*[1]; //Default is 4 faces

    int **inpoel = new int*[1];

    /*Edges of element*/
    //column index is element index, each row in column contents index of edge belong to element, number of row is 4 because of default quad element
    //column index is edge index, each row in column contents index of element which edge is belong to, row 3 contents pointer
    int **inedel = new int*[1];
    int **ineled = new int*[1];

    double **Points = new double*[1]; //array nay resize trong ham loadMesh
    double **normalVector = new double*[1];

    int **Elements1D = new int*[1]; //array nay resize trong ham loadMesh
    int **Elements2D = new int*[1]; //array nay resize trong ham loadMesh
    int **BoundaryType = new int*[1]; //BoundaryType: column 0 is boundary group (from 0), column 1 is boundary type (1, 2, 3, 4), column 3 is boundary method
    int *markPointsAtBC = new int[1];
    int *MasterElemOfEdge = new int[1];

    /*Gauss points on edges*/
    double **edgeGaussPoints_a = new double*[1];
    double **edgeGaussPoints_b = new double*[1];

	/*Vector contents BC edges name and location of them on BC values arrays*/
    //std::vector<int>adressOfBCVals(0, 0);
    int **neighboringElements = new int*[1];

    double **geoCenter = new double*[1];
    double *cellSize = new double[1];
    double *cellArea = new double[1];
    double *localCellSize = new double[1];

    /*derivatives dx/da, dx/db, dy/da, dy/db*/
    //dxa,dxb ... la mang 3D, luu gia tri dao ham dxa, dxb ... tai moi diem Gauss
    double **dxa = new double*[1];
    double **dxb = new double*[1];
    double **dya = new double*[1];
    double **dyb = new double*[1];

    /*Jacobian*/
    //J2D ban dau la mang 3D, luu gia tri Jacobi tai moi diem Gauss
    double **J2D = new double*[1];
    double *J1D;

    /*Array needed for decomposing case*/
    int *Elem2DlocalIdWithRank = new int[1];
    int *rankOf2DElem = new int[1];
    int **PointslocalIdWithRank = new int*[1];

    /*Mesh connection*/
    int **meshConnection = new int*[1];

    /*For Maxwell-Smoluchowsky BC*/
    double **normProjectionOfCenterToBCEdge_realSysCoor = new double*[1];
    double **normProjectionOfCenterToBCEdge_standardSysCoor = new double*[1];
    double *distanceFromCentroidToBCEdge = new double[1];

    double **GaussPtsOnBCEdge_x = new double*[1];
    double **GaussPtsOnBCEdge_y = new double*[1];
    double **distanceFromGaussPtsToCentroid = new double*[1];

    double **GaussPtsOnBCEdge_unitVector_x = new double*[1];
    double **GaussPtsOnBCEdge_unitVector_y = new double*[1];
}

namespace mathVar {
    double *wGaussVol = new double[1];
    double *xGaussVol = new double[1];
    double *wGaussLobattoVol = new double[1];
    double *xGaussLobattoVol = new double[1];

    double *wGaussSur = new double[1];
    double *xGaussSur = new double[1];
    double *wGaussLobattoSur = new double[1];
    double *xGaussLobattoSur = new double[1];

    double *B = new double[1];
    double *dBa = new double[1];
    double *dBb = new double[1];

    //std::vector<std::vector<std::vector<double>>> GaussPts, wGaussPts, GaussLobattoPts, wGaussLobattoPts,
    //    BPts_Quad, dBaPts_Quad, dBbPts_Quad, BPts_Tri, dBaPts_Tri, dBbPts_Tri;
    //Cac mang tren la mang 3D, chua cac gia tri tai tung diem Gauss cua 1 cell
    double **GaussPts = new double*[1];
    double **wGaussPts = new double*[1];
    double **GaussLobattoPts = new double*[1];
    double **wGaussLobattoPts = new double*[1];
    double **BPts_Quad = new double*[1];
    double **dBaPts_Quad = new double*[1];
    double **dBbPts_Quad = new double*[1];
    double **BPts_Tri = new double*[1];
    double **dBaPts_Tri = new double*[1];
    double **dBbPts_Tri = new double*[1];
}

/*Conservative variables declaration*/
double **rho = new double*[1];
double **rhou = new double*[1];
double **rhov = new double*[1];
double **rhoE = new double*[1];

double **rhoN = new double*[1];
double **rhouN = new double*[1];
double **rhovN = new double*[1];
double **rhoEN = new double*[1];

double **rho0 = new double*[1];
double **rhou0 = new double*[1];
double **rhov0 = new double*[1];
double **rhoE0 = new double*[1];

double **rhoResArr = new double*[1];
double **rhouResArr = new double*[1];
double **rhovResArr = new double*[1];
double **rhoEResArr = new double*[1];

/*Primary variables declaration*/
/*
std::vector<std::vector<double>>
u(1, std::vector<double>(1, 0.0)),
v(1, std::vector<double>(1, 0.0)),
e(1, std::vector<double>(1, 0.0)),
p(1, std::vector<double>(1, 0.0)),
T(1, std::vector<double>(1, 0.0)),
mu(1, std::vector<double>(1, 0.0));
*/

/*Auxilary variables*/
namespace BR1Vars {
//X direction
double **rhoX = new double*[1];
double **rhouX = new double*[1];
double **rhovX = new double*[1];
double **rhoEX = new double*[1];

/*Y direction*/
double **rhoY = new double*[1];
double **rhouY = new double*[1];
double **rhovY = new double*[1];
double **rhoEY = new double*[1];
}

namespace BR2Vars {
//X direction
double **rhoXVol = new double*[1];
double **rhouXVol = new double*[1];
double **rhovXVol = new double*[1];
double **rhoEXVol = new double*[1];

/*Y direction*/
double **rhoYVol = new double*[1];
double **rhouYVol = new double*[1];
double **rhovYVol = new double*[1];
double **rhoEYVol = new double*[1];

//X direction
double **rhoXSurMaster = new double*[1];
double **rhouXSurMaster = new double*[1];
double **rhovXSurMaster = new double*[1];
double **rhoEXSurMaster = new double*[1];
double **rhoXSurSlave = new double*[1];
double **rhouXSurSlave = new double*[1];
double **rhovXSurSlave = new double*[1];
double **rhoEXSurSlave = new double*[1];

/*Y direction*/
double **rhoYSurMaster = new double*[1];
double **rhouYSurMaster = new double*[1];
double **rhovYSurMaster = new double*[1];
double **rhoEYSurMaster = new double*[1];
double **rhoYSurSlave = new double*[1];
double **rhouYSurSlave = new double*[1];
double **rhovYSurSlave = new double*[1];
double **rhoEYSurSlave = new double*[1];
}

/*Interface conservative variables*/
namespace surfaceFields {
    double **rho = new double*[1];
    double **rhou = new double*[1];
    double **rhov = new double*[1];
    double **rhoE = new double*[1];

    /*Interface values*/
    //Auxilary equation
    double **aux_rho = new double*[1];
    double **aux_rhou = new double*[1];
    double **aux_rhov = new double*[1];
    double **aux_rhoE = new double*[1];

    //NSF equation
    //inviscid/viscous flux
    double **invis_rho = new double*[1];
    double **invis_rhou = new double*[1];
    double **invis_rhov = new double*[1];
    double **invis_rhoE = new double*[1];

    double **Vis_rho = new double*[1];
    double **Vis_rhou = new double*[1];
    double **Vis_rhov = new double*[1];
    double **Vis_rhoE = new double*[1];

    double **T = new double*[1];

    //Derivatives
    //Ox
    double **dRhoX = new double*[1];
    double **dRhouX = new double*[1];
    double **dRhovX = new double*[1];
    double **dRhoEX = new double*[1];
    //Oy
    double **dRhoY = new double*[1];
    double **dRhouY = new double*[1];
    double **dRhovY = new double*[1];
    double **dRhoEY = new double*[1];
}

namespace volumeFields {
    //Volume values
    //Cac array cua volumeFields la 3D array
    double **rhoVolGauss = new double*[1];
    double **rhouVolGauss = new double*[1];
    double **rhovVolGauss = new double*[1];
    double **rhoEVolGauss = new double*[1];
    double **drhoXVolGauss = new double*[1];
    double **drhoYVolGauss = new double*[1];
    double **T = new double*[1];
}

//Lax-Friedrich constant
double *LxFConst = new double[1];
double *DiffusiveFluxConst = new double[1];

//StiffMatrixCoefficients
double **stiffMatrixCoeffs = new double*[1];

//Limiting coefficients
double *theta1Arr = new double[1];
double *theta2Arr = new double[1];

namespace SurfaceBCFields
{
    //int **BCPointsInfor;
    std::vector<std::vector<int>> BCPointsInfor;

    //double *TBc = new double[1];
    //double *uBc = new double[1];
    //double *vBc = new double[1];
    //double *pBc = new double[1];

    double **TBc = new double*[1];
    double **uBc = new double*[1];
    double **vBc = new double*[1];
    double **pBc = new double*[1];
    double **rhoBc = new double*[1];

    int *localGlobalBCEdgesMatching = new int[1];

    /*
    //Derivatives
    //Ox
    double **GaussDRhoX = new double*[1];
    double **GaussDRhouX = new double*[1];
    double **GaussDRhovX = new double*[1];
    double **GaussDRhoEX = new double*[1];
    //Oy
    double **GaussDRhoY = new double*[1];
    double **GaussDRhouY = new double*[1];
    double **GaussDRhovY = new double*[1];
    double **GaussDRhoEY = new double*[1];
    */
}

namespace limitVal {
    bool *troubleCellsMarker = new bool[1];
}

//for debugging
namespace debug
{
	//std::vector<double>
		//minRhoArr(1, 1.0),
		//minRhoeArr(1, 1.0);
}
