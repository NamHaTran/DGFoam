#ifndef SUPPORTFUNCS_WAVETRANSIMISSIVE_H
#define SUPPORTFUNCS_WAVETRANSIMISSIVE_H
#include <vector>
#include <sstream>
namespace waveTransmissive
{
    extern bool *includeSoundSpeed_p;
    extern bool *includeSoundSpeed_T;
    extern bool *includeSoundSpeed_u;

    void readCondition(int bcGrp, std::ifstream &FileFlux, std::string file);
}
#endif // SUPPORTFUNCS_WAVETRANSIMISSIVE_H
