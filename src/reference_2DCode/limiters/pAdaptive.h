#ifndef PADAPTIVE_H
#define PADAPTIVE_H
#include <vector>

namespace limiter
{
    namespace pAdaptive
    {
        void limiter();

        void limiter_1Elem(int element, int valType);

        std::vector<double> pAdaptiveChildFunction_Quad(int element, int valType, double IPlus, double IMinus, double JPlus, double JMinus);

        std::vector<double> pAdaptiveChildFunction_Tri(int element, int valType, double IPlus, double IMinus, double JMinus);
    }
}
#endif // PADAPTIVE_H
