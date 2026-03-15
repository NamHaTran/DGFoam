#ifndef NONEQMBCSGENFUNCS_H
#define NONEQMBCSGENFUNCS_H
#include <string>

namespace nonEquilibriumBCs {

    /**
     * @brief Function reads surface data (T, u, v) in file TSurface.txt, uSurface.txt and vSurface.txt.
     * @param Loc: location of file
     */
    void readSurfaceValues(std::string Loc);

    /**
     * @brief Function writes surface data (T, u, v) to file TSurface.txt, uSurface.txt and vSurface.txt.
     * @param Loc: location of file
     */
    void writeSurfaceValues(std::string Loc);

    /**
     * @brief Function runs all time varying BCs and update results to SurfaceBCFields arrays.
     *
     * Functions runs after TVD_RK3 function finished.
     */
    void updateBCs();
}
#endif // NONEQMBCSGENFUNCS_H
