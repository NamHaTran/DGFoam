#include "limiterController.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "./limiters/positivityPreserving/positivityPreserving.h"
#include "./limiters/massDiffusion/massDiffusion.h"
#include "./limiters/pAdaptive/pAdaptive.h"

#include <sstream>
#include <fstream>
#include <string>

namespace limiter
{
    //ham limiter chay sau khi hoan thanh 1 step TVDRK3
    void limiter_1InnerStep()
    {
        if (mathVar::orderElem != 0)
        {
            limiter::pAdaptive::limiter();

            limiter::positivityPreserving::limiter();
        }
    }

    void limiter_1OutterStep()
    {
        if (limitVal::massDiffusion)
        {
            limiter::massDiffusion::limiter();
        }
    }

    /**
     * @brief Function limits Rho by applying Positivity Preserving Limiter.
     */
    void limitRho_PositivityPreserving()
    {
        if (mathVar::orderElem!=0)
        {
            if (limitVal::PositivityPreserving)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {
                    theta1Arr[nelem] = limiter::positivityPreserving::calcTheta1(nelem);
                }
            }
        }
    }

    namespace IOLimiter {
        void readSelectedLimiters()
        {
            //Read selected limiters
            std::string FileDir(systemVar::wD + "/CASES/" + systemVar::caseName + "/System"), fileName("LimiterSettings.txt");
            std::string FileLoc(FileDir + "/" + fileName);
            std::ifstream FileFlux(FileLoc.c_str());
            if (FileFlux)
            {
                std::string line, Word;
                while (std::getline(FileFlux, line))
                {
                    std::istringstream line2str(line);
                    std::vector<std::string> ptr;
                    //Split <line2str> stringstream into array of words
                    while ((line2str >> Word))
                    {
                        ptr.push_back(Word);
                    }

                    int numWrd = static_cast<int>(ptr.size());
                    if (numWrd >= 2)
                    {
                        std::istringstream strdata(ptr[1]);

                        if (ptr[0].compare("limiter") == 0) //get selected limiter(s)
                        {
                            for (int i = 0; i < numWrd - 1; i++)
                            {
                                limitVal::limiterName.push_back(ptr[i + 1]);
                            }

                            if (limitVal::limiterName.size() > 0)
                            {
                                for (int ilimiter = 0; ilimiter < static_cast<int>(limitVal::limiterName.size()); ilimiter++)
                                {
                                    if (limitVal::limiterName[ilimiter].compare("positivityPreserving") == 0)
                                    {
                                        limitVal::PositivityPreserving = true;
                                    }
                                    if (limitVal::limiterName[ilimiter].compare("pAdaptive") == 0)
                                    {
                                        limitVal::PAdaptive = true;
                                    }
                                    if (limitVal::limiterName[ilimiter].compare("massDiffusion") == 0)
                                    {
                                        limitVal::massDiffusion = true;
                                    }
                                }
                            }

                            break;
                        }
                    }
                }
            }
            else
            {
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
            }

            //Read settings of each limiter
			if (limitVal::massDiffusion)
				limiter::IOMassDiff::readSetting();
				
			if (limitVal::PositivityPreserving)
				limiter::IOPositivity::readSetting();
        }
    }
}
