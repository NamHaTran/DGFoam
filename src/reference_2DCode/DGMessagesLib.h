#ifndef DGMESSAGESLIB_H_INCLUDED
#define DGMESSAGESLIB_H_INCLUDED
#include <string>
namespace message
{
	/*Function create header of DG's data files*/
	std::string headerFile();

	/*Function display undefined keyWord error*/
	std::string undfKeyW(std::string keyW, std::string location);

	/*Function display opening file error*/
	std::string opFError(std::string fileName, std::string location);

	/*Function writes logFile to report error to user*/
	void writeLog(std::string location, std::string caseName, std::string str);

	/*Functions gets time data from system*/
	std::string getTime();

	/*Function display undefined boundary condition type error*/
	std::string undfBcType(std::string keyW, std::string fileName, std::string bcType);

	/*Function create header of files p T U*/
	std::string headerpTU(std::string file);

	/*Function display error of Slip Bc uncompatibility*/
	std::string SlipBcCompatibleError(int bcGrp);

	/*Function display error of Bc uncompatibility*/
	std::string BcCompatibleError(int edgeGrp);

	/*Function display error of wrong values of orderElem and nGauss*/
	std::string nGaussOrderElemError();

	/*Help functions--------------------------------*/

	/*Function display help for unvReader ultility*/
	void UnvReaderHelp();

	/*Function display help for boundary conditions*/
	void BCsHelp();

    /*Function prints out case's informations and asks user whether to continue or not*/
    void showCaseInformations();

    void showWarning();
}

/*Function displays error message and exit program*/
void exitDG(std::string str);
#endif // DGMESSAGESLIB_H_INCLUDED
