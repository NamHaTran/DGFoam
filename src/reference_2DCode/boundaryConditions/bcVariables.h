 #ifndef BCVARIABLES_H
#define BCVARIABLES_H
#include <vector>
namespace BCVars {

    extern std::vector<bool> DirichletAppMethUStrong, NewmannAppMethGradUStrong,
    DirichletAppMethTStrong, NewmannAppMethGradTStrong,
    DirichletAppMethPStrong, NewmannAppMethGradPStrong,
    DirichletAppMethGeneralBCStrong, NewmannAppMethGradGeneralBCStrong;

    namespace generalBCId {
        static constexpr int symmetry = 7;
        static constexpr int matched = 10;
    }

    namespace pressureBCId {
        static constexpr int interpFrmDensity = 11;
        static constexpr int zeroGrad = 2;
        static constexpr int fixedValue = 1;
        static constexpr int inletOutlet = 4;
        static constexpr int zeroGradRhoUncorrectP = 12;

        //Custom BC reflectGradRho
        static constexpr int reflectGradRho = 13;

        //Custom BC interiorSide (add 22/07/2021)
        static constexpr int interiorSide = 14;

        //Custom BC reflectGradRho (add 23/09/2021)
        static constexpr int zeroRhoGrad = 15;

        //Custom BC reflectGradRho (add 03/06/2022)
        static constexpr int waveTransmissive = 16;
    }

    namespace temperatureBCId {
        static constexpr int temperatureJump = 11;
        static constexpr int zeroGrad = 2;
        static constexpr int fixedValue = 1;
        static constexpr int inletOutlet = 4;

        //Custom BC interiorSide (add 22/07/2021)
        static constexpr int interiorSide = 14;

        //Custom BC SmoluchowskyTJump (add 22/07/2021)
        static constexpr int SmoluchowskyTJump = 15;

        //Custom BC reflectGradRho (add 03/06/2022)
        static constexpr int waveTransmissive = 16;
    }

    namespace velocityBCId {
        static constexpr int fixedValue = 1;
        static constexpr int noSlip = 2;
        static constexpr int movingWall = 3;
        static constexpr int slipWall = 5;
        static constexpr int inletOutlet = 4;
        static constexpr int zeroGrad = 6;

        //Custom BC interiorSide (add 22/07/2021)
        static constexpr int interiorSide = 14;

        //Custom BC MaxwellSlip (add 19/07/2021)
        static constexpr int MaxwellSlip = 15;

        //Custom BC reflectGradRho (add 03/06/2022)
        static constexpr int waveTransmissive = 16;
    }
}
#endif // BCVARIABLES_H
