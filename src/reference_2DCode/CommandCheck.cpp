#include "CommandCheck.h"
#include <string>

/*! \brief postProcessKey.
 *         Namespace contains function of checking post-processing commands.
 *
 *  No detailed description.
 */
namespace postProcessKey
{
	/*Function return true if Exit is available*/
	bool checkExit(std::string cmd)
	{
		bool trigger(false);
        if ((cmd.compare("Exit") == 0) || (cmd.compare("exit") == 0) || (cmd.compare("EXIT") == 0) || (cmd.compare("End") == 0) || (cmd.compare("end") == 0) || (cmd.compare("END") == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}

/*! \brief preProcessKey.
 *         Namespace contains function of checking pre-processing commands.
 *
 *  No detailed description.
 */
namespace preProcessKey
{
    /**
     * @brief Function return true if UnvToDG is available.
     *
     * UnvToDG is used for converting Unv mesh to DG mesh.
     *
     * @param cmd command inputed from terminal.
     * @return trigger
     */
	bool checkUnvReader(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("UnvToDG") == 0) || (cmd.compare("unvtodg") == 0) || (cmd.compare("UNVTODG") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

    /**
     * @brief Function return true if UnvReaderHelp is available.
     *
     * UnvReaderHelp shows helps of UnvReader function.
     *
     * @param cmd command inputed from terminal.
     * @return trigger
     */
	bool checkUnvHelper(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("UnvReaderHelp") == 0) || (cmd.compare("unvreaderhelp") == 0) || (cmd.compare("unvreaderhelp") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

    /**
     * @brief Function return true if bcHelp is available.
     *
     * bcHelp shows helps of available boundary conditions.
     *
     * @param cmd command inputed from terminal.
     * @return trigger
     */
	bool checkBCsHelper(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("bchelp") == 0) || (cmd.compare("BcHelp") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

    /**
     * @brief Function return true if reSubmit is available.
     *
     * reSubmit reread file submitCase.txt to submit new case.
     *
     * @param cmd
     * @return trigger
     */
	bool reSubmit(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("resubmit") == 0) || (cmd.compare("reSubmit") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

    /**
     * @brief Function return true if mappResults is available.
     *
     * mappResults maps results from low order case to higher order case.
     *
     * @param cmd
     * @return trigger
     */
	bool mappResults(std::string cmd)
	{
		bool trigger(false);
        if ((cmd.compare("mapresults") == 0) || (cmd.compare("mapResults") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

    /**
     * @brief Function return true if exportMeshToMetis is available.
     *
     * exportMeshToMetis exports DG mesh to Metis's format. Metis is a tool to distribute mesh into to defined number of cores.
     *
     * @param cmd
     * @return trigger
     */
    bool exportMeshToMetis(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("DGMeshToMetis") == 0) || (cmd.compare("dgmeshtometis") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    /**
     * @brief Function return true if testMeshPartitionResult is available.
     *
     * testMeshPartitionResult exports partitioned mesh to techplot format which can be opened on Paraview.
     *
     * @param cmd
     * @return trigger
     */
    bool testMeshPartitionResult(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("testmeshpartitionresult") == 0) || (cmd.compare("testMeshPartitionResult") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    /**
     * @brief Function return true if decomposeCase is available.
     *
     * decomposeCase decomposes mesh and boundary condtions of case to all cores.
     *
     * @param cmd
     * @return trigger
     */
    bool decomposeCase(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("decomposecase") == 0) || (cmd.compare("decomposeCase") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    /**
     * @brief Function return true if reconstructLatestTime is available.
     *
     * reconstructLatestTime reconstruct mesh and results of decomposed case.
     *
     * @param cmd
     * @return trigger
     */
    bool reconstructLatestTime(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("reconstructCase") == 0) || (cmd.compare("reconstructcase") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    /**
     * @brief Function return true if checkPartitionedMesh is available.
     *
     * checkPartitionedMesh loops over all mesh cells to find isolated cells, then changes the master core to the lowest index master core of neighbour cell.\n
     * Explainations:
     * - Isolated cell is the cell which has ALL neighbour cells belong to different master cores in comparison with itself.
     * - Master core of cell is the core containing the cell after decomposing mesh.
     *
     * @param cmd
     * @return trigger
     */
    bool checkPartitionedMesh(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("checkParMesh") == 0) || (cmd.compare("checkparmesh") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

	namespace debug
	{
		bool checkElement(std::string cmd)
		{
			bool trigger(false);
			if ((cmd.compare("checkElementInfor") == 0) || (cmd.compare("checkelementinfor") == 0))
			{
				trigger = true;
			}
			return trigger;
		}
	}
}

/*! \brief processKey.
 *         Namespace contains function of checking processing commands.
 *
 *  No detailed description.
 */
namespace processKey
{
    /**
     * @brief Function return true if checkDGRun is available.
     *
     * @param cmd
     * @return trigger
     */
	bool checkDGRun(std::string cmd)
	{
		bool trigger(false);
        if ((cmd.compare("DG2D") == 0) || (cmd.compare("dg2d") == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}
