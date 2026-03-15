#include "nonEqmBCs_Vars.h"

namespace nonEqmSurfaceField {
    double **TBc = new double*[1];
    double **uBc = new double*[1];
    double **vBc = new double*[1];
    double **pBc = new double*[1];
}

namespace meshVar {
    double **GaussPtsOnBCEdge_x = new double*[1];
    double **GaussPtsOnBCEdge_y = new double*[1];
    double **distanceFromGaussPtsToCentroid = new double*[1];

    double **GaussPtsOnBCEdge_unitVector_x = new double*[1];
    double **GaussPtsOnBCEdge_unitVector_y = new double*[1];
}
