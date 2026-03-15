#ifndef COMMANDCHECK_H_INCLUDED
#define COMMANDCHECK_H_INCLUDED
#include <string>

namespace postProcessKey
{
	/*Function return true if Exit is available*/
	bool checkExit(std::string cmd);
}

namespace preProcessKey
{
	/*Function return true if UnvToDG is available*/
	bool checkUnvReader(std::string cmd);

	/*Function return true if UnvReaderHelp is available*/
	bool checkUnvHelper(std::string cmd);

	/*Function return true if bcHelp is available*/
	bool checkBCsHelper(std::string cmd);

	/*Function return true if reSubmit is available*/
	bool reSubmit(std::string cmd);

	bool mappResults(std::string cmd);

    bool exportMeshToMetis(std::string cmd);

    bool testMeshPartitionResult(std::string cmd);

    bool decomposeCase(std::string cmd);

    bool reconstructLatestTime(std::string cmd);

    bool checkPartitionedMesh(std::string cmd);

	namespace debug
	{
		bool checkElement(std::string cmd);
	}
}

namespace processKey
{
	/*Function return true if DG2D is available*/
	bool checkDGRun(std::string cmd);
}

#endif // COMMANDCHECK_H_INCLUDED
