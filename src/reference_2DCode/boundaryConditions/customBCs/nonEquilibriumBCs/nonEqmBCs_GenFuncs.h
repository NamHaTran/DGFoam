#ifndef NONEQMBCS_GENFUNCS_H
#define NONEQMBCS_GENFUNCS_H
#include <string>

namespace nonEquilibriumBCs {
    /**
     * @brief Function runs all time varying BCs and update results to SurfaceBCFields arrays.
     *
     * Functions runs after TVD_RK3 function finished.
     */
    void updateBCs();

    void resizeSurfaceFields();
}

namespace nonEquilibriumBCs_IO {
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

    void writeToFile(std::string loc, std::string name, double **array);

    void readFromFile(std::string loc, std::string fileName, bool exitWhenFileNotFound, double **array);

    void writeGaussPtsCoor(std::string loc, std::string name);
}
#endif // NONEQMBCS_GENFUNCS_H
