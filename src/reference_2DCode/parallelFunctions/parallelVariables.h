#ifndef PARALLELVARIABLES_H
#define PARALLELVARIABLES_H
namespace parallelBuffer {
    extern double **xCoor;
    extern double **yCoor;
    //extern double **aCoor;
    //extern double **bCoor;

    //extern double *theta1;
    //extern double *theta2;
    extern int *elemType;

    namespace surfaceGaussPt {
        //extern double **rho;
        //extern double **rhou;
        //extern double **rhov;
        //extern double **rhoE;

    /*
        extern double **drhoX;
        extern double **drhouX;
        extern double **drhovX;
        extern double **drhoEX;
        extern double **drhoY;
        extern double **drhouY;
        extern double **drhovY;
        extern double **drhoEY;*/

        //extern double **T;
    }
}
#endif // PARALLELVARIABLES_H
