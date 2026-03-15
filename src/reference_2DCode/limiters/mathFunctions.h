#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H
#include <vector>

namespace limiter
{
    //MATHEMATIC FUNCTIONS FOR LIMITER
    namespace mathForLimiter
    {
        double minmod(std::vector<double> vectorArguments);

        double modifiedMinmod(std::vector<double> vectorArguments, double condition);

        int getSignOfDouble(double input);

        double calM(int element, int valType);

        void getNeighborElements();
    }
}
#endif // MATHFUNCTIONS_H
