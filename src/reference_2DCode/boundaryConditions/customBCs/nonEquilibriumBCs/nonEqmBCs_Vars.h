#ifndef NONEQMBCS_VARS_H
#define NONEQMBCS_VARS_H

namespace nonEqmSurfaceField {
    extern double **TBc;
    extern double **uBc;
    extern double **vBc;
    extern double **pBc;
}

namespace meshVar {
    extern double **GaussPtsOnBCEdge_x;
    extern double **GaussPtsOnBCEdge_y;
    extern double **distanceFromGaussPtsToCentroid;

    extern double **GaussPtsOnBCEdge_unitVector_x;
    extern double **GaussPtsOnBCEdge_unitVector_y;
}
#endif // NONEQMBCS_VARS_H
