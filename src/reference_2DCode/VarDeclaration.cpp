#include "VarDeclaration.h"
#include <string>
#include "ConstDeclaration.h"
#include <vector>
#include "DGMessagesLib.h"

namespace systemVar
{
	std::string wD;
	std::string caseName;
	std::string pwd;
	std::string cmd;
	bool endKey(false);

    double CFL(0.5); //Courant number
	double Ttime(0.0); //Total time
    int wrtI(0); //write interval
    bool wrtLog(true), loadSavedCase(false);

    int ddtScheme(1);
    double epsilon(1e-13);

    int iterCount(0), savingCout(0), loadConstCount(0);
	double rhoResNorm(1.0), rhouResNorm(1.0), rhovResNorm(1.0), rhoEResNorm(1.0);

	bool initializedOrNot(false), runPreProcess(false);

    //1: BR1, 2: BR2
    int auxVariables(1);

    //For parallel computing
    int totalProc(4), currentProc(0);
    bool runDecomposeCaseFnc(false), parallelMode(true), runCheckParMesh(false);

    //Detect first iteration
    bool firstIter(true);
	
	std::string headerFile(message::headerFile());

    //sendRecvOrder
    int **sendRecvOrder = new int*[1];
    int sendRecvOrder_length;

    //Read/write mode
    std::string readWriteMode(" ");
}

namespace meshVar
{
    int nelem1D(0), nelem2D(0), npoin(0), nBc(0);

	/*Default value*/
	//number of nodes per element
    int const nnode(4);
	//number of edges per element (default is 4)
	int const nedel(4);  //change nedel value at file .h

    //Cac array o day dung kieu std::vector vi size thay doi trong qua trinh doc luoi
	/*Elements surrounding point*/
	std::vector<int> esup1, esup2;

	/*Points surrounding point*/
	std::vector<int> psup1, psup2;

	/*Variables help to save mesh data*/
	int inpoedCount(0);  //can be used for normalVector, MasterElemOfEdge, ineled

	int numBCEdges(0);

    //Number of BC groups
    int numBCGrp(0);
}

namespace mathVar
{
    int orderElem(0), orderOfAccuracy(0);
    bool solveTFailed(false);

    int nGauss1D(0), nGauss2D(0);
}

namespace material
{
    double gamma(1.4), R(287.0), Pr(0.72), Cp(0.0), Cv(0.0);

    namespace viscosityCoeff {
        namespace Sutherland {
            double
            As(0.001),
            Ts(110.4);
        }

        namespace powerLaw_VHS {
            double
            molMass(0.0),
            kBoltzmann(1.380649e-23),
            omega(0.0),
            TRef(0.0),
            dRef(3.595e-10);
        }

        namespace constant {
            double mu(0.0);
        }
    }

    namespace viscousityModel {
        bool constant(false),
        power_VHS(false),
        sutherland(false);
    }
}

namespace iniValues
{
	double uIni(0.0), vIni(0.0), wIni(0.0), pIni(0.0), TIni(0.0), muIni(0.0), rhoIni(0.0), eIni(0.0);
}

namespace bcValues
{
    double uBCFixed[bcSize] = { 0.0 },
        vBCFixed[bcSize] = { 0.0 },
        wBCFixed[bcSize] = { 0.0 },
        pBCFixed[bcSize] = { 0.0 },
        TBCFixed[bcSize] = { 0.0 };

	/*Variables help identify boundary condition at each group*/
	int UBcType[bcSize] = {};
	int pBcType[bcSize] = {};
	int TBcType[bcSize] = {};

    /*Variables of slip & temperature jump boundary conditions*/
    double sigmaU(1.0), sigmaT(1.0);

    /*Flags of time varying non equilibrium BCs*/
    //U
    bool slipBCFlag(false);
    //T
    bool temperatureJump(false);

    /*Flag of time varying BCs*/
    bool timeVaryingBCAvailable(false);
}

namespace flowProperties
{
    bool viscous(true), massDiffusion(false);
    double Mach(0.0);
    bool subsonic(true);
}

//time step
double dt(1e-20);
double runTime(0.0);

namespace limitVal
{
    double TUp(5000), TDwn(20);
	double rhoUp(12.5), rhoDwn(0.0001);
	double rhoEUp(0.0), rhoEDwn(0.0);  //computed at initialValue function
	bool limitTOrNot(false), limitFlagLocal(false), limitFlagGlobal(false);
	int numOfLimitCell(0);

	namespace pAdaptive
	{
		bool limitFlagLocal(false), limitFlagGlobal(false);
		int numOfLimitCell(0);
	}

	std::vector<std::string> limiterName;

    //keys of limiter
    bool PositivityPreserving(false), PAdaptive(false), massDiffusion(false),
        runningMassDiffLimiter(false) //bien kiem tra xem mass diffusion limiter co dang chay hay khong
    ;
	namespace PositivityPreservingSettings
	{
		/*version:
		- 1: full
		- 2: simplified
		*/
		int version(2);
	}
}

namespace controlFlag {
    namespace sequence {
    bool checkUnvReader(false),
    checkBCsHelper(false),
    checkUnvHelper(false),
    reSubmit(false),
    exportMeshToMetis(false),
    testMeshPartitionResult(false),
    debug_checkElement(false),
    decomposeCase(false),
    reconstructLatestTime(false),
    checkPartitionedMesh(false);
    }

    namespace parallel {
    bool checkDGRun(false),
    mapResults(false);
    }

    /* Cac flag trong namespace nay check xem folder chua kq co file TSurface, uSurface ... hay khong
     * Neu khong co, cac gia tri cua cac field nay set ve gia tri mac dinh o dieu kien bien
    */
    namespace fileAvailFlags {
    bool fileTSurface(false),
        fileuSurface(false),
        filevSurface(false);
    }
}

namespace warningFlag {
    bool reversedFlowOccur(false);
}

namespace DGSchemes {
    //Dung cho truong hop mass diffusion, true neu giai T implicit, false neu giai T explicit
    bool solveTImplicitly(true);

    namespace fluxControl {
        //Tai version hien tai, flux cua bien phu va viscous flux la central flux, flux control chi apply cho
        //convective flux
        bool LxF(false), Roe(false), HLL(false), HLLC(false), central(false);
    }
}

namespace debugVars {
int element(0);
}
