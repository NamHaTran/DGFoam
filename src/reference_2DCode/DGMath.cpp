#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMessagesLib.h"
#include <math.h>
#include <tuple>  //Include this for returning multiple values in function
#include <algorithm>
#include <iostream>

#include "./parallelFunctions/parallelVariables.h"

#include "debuggingFuncs.h"

namespace math
{
void getGaussCoeffs(double* xGauss, double* wGauss, int nGauss)
    {
        if (nGauss==0)
        {
            xGauss[nGauss] = 0.0;
            wGauss[nGauss] = 2.0;
        }
        else if (nGauss==1)
        {
            xGauss[nGauss - 1] = -0.577350269189625764509148780502;
            xGauss[nGauss] = 0.577350269189625764509148780502;

            wGauss[nGauss - 1] = 1.0;
            wGauss[nGauss] = 1.0;
        }
        else if (nGauss==2)
        {
            xGauss[nGauss - 2] = -0.774596669241483377035853079956;
            xGauss[nGauss - 1] = 0.0;
            xGauss[nGauss] = 0.774596669241483377035853079956;

            wGauss[nGauss - 2] = 0.555555555555555555555;
            wGauss[nGauss - 1] = 0.888888888888888888888;
            wGauss[nGauss] = 0.55555555555555555555555;
        }
        else if (nGauss==3)
        {
            xGauss[nGauss - 3] = -0.8611363115940526;
            xGauss[nGauss - 2] = -0.3399810435848563;
            xGauss[nGauss - 1] = 0.3399810435848563;
            xGauss[nGauss] = 0.8611363115940526;

            wGauss[nGauss - 3] = 0.3478548451374538;
            wGauss[nGauss - 2] = 0.6521451548625461;
            wGauss[nGauss - 1] = 0.6521451548625461;
            wGauss[nGauss] = 0.3478548451374538;
        }
        else if (nGauss==4)
        {
            xGauss[nGauss - 4] = -0.9061798459386640;
            xGauss[nGauss - 3] = -0.5384693101056831;
            xGauss[nGauss - 2] = 0.0000000000000000;
            xGauss[nGauss - 1] = 0.5384693101056831;
            xGauss[nGauss] = 0.9061798459386640;

            wGauss[nGauss - 4] = 0.2369268850561891;
            wGauss[nGauss - 3] = 0.4786286704993665;
            wGauss[nGauss - 2] = 0.5688888888888889;
            wGauss[nGauss - 1] = 0.4786286704993665;
            wGauss[nGauss] = 0.2369268850561891;
        }
        else if (nGauss==5)
        {
            xGauss[nGauss - 5] = -0.9324695142031521;
            xGauss[nGauss - 4] = -0.6612093864662645;
            xGauss[nGauss - 3] = -0.2386191860831969;
            xGauss[nGauss - 2] = 0.2386191860831969;
            xGauss[nGauss - 1] = 0.6612093864662645;
            xGauss[nGauss] = 0.9324695142031521;

            wGauss[nGauss - 5] = 0.1713244923791704;
            wGauss[nGauss - 4] = 0.3607615730481386;
            wGauss[nGauss - 3] = 0.4679139345726910;
            wGauss[nGauss - 2] = 0.4679139345726910;
            wGauss[nGauss - 1] = 0.3607615730481386;
            wGauss[nGauss] = 0.1713244923791704;
        }
        else if (nGauss==6)
        {
            xGauss[nGauss - 6] = -0.9491079123427585;
            xGauss[nGauss - 5] = -0.7415311855993945;
            xGauss[nGauss - 4] = -0.4058451513773972;
            xGauss[nGauss - 3] = 0.0;
            xGauss[nGauss - 2] = 0.4058451513773972;
            xGauss[nGauss - 1] = 0.7415311855993945;
            xGauss[nGauss] = 0.9491079123427585;

            wGauss[nGauss - 6] = 0.1294849661688697;
            wGauss[nGauss - 5] = 0.2797053914892766;
            wGauss[nGauss - 4] = 0.3818300505051189;
            wGauss[nGauss - 3] = 0.4179591836734694;
            wGauss[nGauss - 2] = 0.3818300505051189;
            wGauss[nGauss - 1] = 0.2797053914892766;
            wGauss[nGauss] = 0.1294849661688697;
        }
        else if (nGauss==7)
        {
            xGauss[nGauss - 7] = -0.9602898564975363;
            xGauss[nGauss - 6] = -0.7966664774136267;
            xGauss[nGauss - 5] = -0.5255324099163290;
            xGauss[nGauss - 4] = -0.1834346424956498;
            xGauss[nGauss - 3] = 0.1834346424956498;
            xGauss[nGauss - 2] = 0.5255324099163290;
            xGauss[nGauss - 1] = 0.7966664774136267;
            xGauss[nGauss] = 0.9602898564975363;

            wGauss[nGauss - 7] = 0.1012285362903763;
            wGauss[nGauss - 6] = 0.2223810344533745;
            wGauss[nGauss - 5] = 0.3137066458778873;
            wGauss[nGauss - 4] = 0.3626837833783620;
            wGauss[nGauss - 3] = 0.3626837833783620;
            wGauss[nGauss - 2] = 0.3137066458778873;
            wGauss[nGauss - 1] = 0.2223810344533745;
            wGauss[nGauss] = 0.1012285362903763;
        }
        else if (nGauss==8)
        {
            xGauss[nGauss - 8] = -0.9681602395076261;
            xGauss[nGauss - 7] = -0.8360311073266358;
            xGauss[nGauss - 6] = -0.6133714327005904;
            xGauss[nGauss - 5] = -0.3242534234038089;
            xGauss[nGauss - 4] = 0.0;
            xGauss[nGauss - 3] = 0.3242534234038089;
            xGauss[nGauss - 2] = 0.6133714327005904;
            xGauss[nGauss - 1] = 0.8360311073266358;
            xGauss[nGauss] = 0.9681602395076261;

            wGauss[nGauss - 8] = 0.0812743883615744;
            wGauss[nGauss - 7] = 0.1806481606948574;
            wGauss[nGauss - 6] = 0.2606106964029354;
            wGauss[nGauss - 5] = 0.3123470770400029;
            wGauss[nGauss - 4] = 0.3302393550012598;
            wGauss[nGauss - 3] = 0.3123470770400029;
            wGauss[nGauss - 2] = 0.2606106964029354;
            wGauss[nGauss - 1] = 0.1806481606948574;
            wGauss[nGauss] = 0.0812743883615744;
        }
        else if (nGauss==9)
        {
            xGauss[nGauss - 9] = -0.9739065285171717;
            xGauss[nGauss - 8] = -0.8650633666889845;
            xGauss[nGauss - 7] = -0.6794095682990244;
            xGauss[nGauss - 6] = -0.4333953941292472;
            xGauss[nGauss - 5] = -0.1488743389816312;
            xGauss[nGauss - 4] = 0.1488743389816312;
            xGauss[nGauss - 3] = 0.4333953941292472;
            xGauss[nGauss - 2] = 0.6794095682990244;
            xGauss[nGauss - 1] = 0.8650633666889845;
            xGauss[nGauss] = 0.9739065285171717;

            wGauss[nGauss - 9] = 0.0666713443086881;
            wGauss[nGauss - 8] = 0.1494513491505806;
            wGauss[nGauss - 7] = 0.2190863625159820;
            wGauss[nGauss - 6] = 0.2692667193099963;
            wGauss[nGauss - 5] = 0.2955242247147529;
            wGauss[nGauss - 4] = 0.2955242247147529;
            wGauss[nGauss - 3] = 0.2692667193099963;
            wGauss[nGauss - 2] = 0.2190863625159820;
            wGauss[nGauss - 1] = 0.1494513491505806;
            wGauss[nGauss] = 0.0666713443086881;
        }
        else if (nGauss==10)
        {
            xGauss[nGauss - 10] = -0.9782286581460570;
            xGauss[nGauss - 9] = -0.8870625997680953;
            xGauss[nGauss - 8] = -0.7301520055740494;
            xGauss[nGauss - 7] = -0.5190961292068118;
            xGauss[nGauss - 6] = -0.2695431559523450;
            xGauss[nGauss - 5] = 0.0;
            xGauss[nGauss - 4] = 0.2695431559523450;
            xGauss[nGauss - 3] = 0.5190961292068118;
            xGauss[nGauss - 2] = 0.7301520055740494;
            xGauss[nGauss - 1] = 0.8870625997680953;
            xGauss[nGauss] = 0.9782286581460570;

            wGauss[nGauss - 10] = 0.0556685671161737;
            wGauss[nGauss - 9] = 0.1255803694649046;
            wGauss[nGauss - 8] = 0.1862902109277343;
            wGauss[nGauss - 7] = 0.2331937645919905;
            wGauss[nGauss - 6] = 0.2628045445102467;
            wGauss[nGauss - 5] = 0.2729250867779006;
            wGauss[nGauss - 4] = 0.2628045445102467;
            wGauss[nGauss - 3] = 0.2331937645919905;
            wGauss[nGauss - 2] = 0.1862902109277343;
            wGauss[nGauss - 1] = 0.1255803694649046;
            wGauss[nGauss] = 0.0556685671161737;
        }
        else if (nGauss==11)
        {
            xGauss[nGauss - 11] = -0.9815606342467192;
            xGauss[nGauss - 10] = -0.9041172563704749;
            xGauss[nGauss - 9] = -0.7699026741943047;
            xGauss[nGauss - 8] = -0.5873179542866175;
            xGauss[nGauss - 7] = -0.3678314989981802;
            xGauss[nGauss - 6] = -0.1252334085114689;
            xGauss[nGauss - 5] = 0.1252334085114689;
            xGauss[nGauss - 4] = 0.3678314989981802;
            xGauss[nGauss - 3] = 0.5873179542866175;
            xGauss[nGauss - 2] = 0.7699026741943047;
            xGauss[nGauss - 1] = 0.9041172563704749;
            xGauss[nGauss] = 0.9815606342467192;

            wGauss[nGauss - 11] = 0.0471753363865118;
            wGauss[nGauss - 10] = 0.1069393259953184;
            wGauss[nGauss - 9] = 0.1600783285433462;
            wGauss[nGauss - 8] = 0.2031674267230659;
            wGauss[nGauss - 7] = 0.2331937645919905;
            wGauss[nGauss - 6] = 0.2491470458134028;
            wGauss[nGauss - 5] = 0.2491470458134028;
            wGauss[nGauss - 4] = 0.2334925365383548;
            wGauss[nGauss - 3] = 0.2031674267230659;
            wGauss[nGauss - 2] = 0.1600783285433462;
            wGauss[nGauss - 1] = 0.1069393259953184;
            wGauss[nGauss] = 0.0471753363865118;
        }
    }

    void getGaussLobattoCoeffs(double* xGaussLobatto, double* wGaussLobatto, int nGauss)
    {
        if (nGauss == 0)
        {
            xGaussLobatto[nGauss] = 0.0;
            wGaussLobatto[nGauss] = 2.0;
        }
        else if (nGauss == 1)
        {
            xGaussLobatto[nGauss - 1] = -1.0;
            xGaussLobatto[nGauss] = 1.0;

            wGaussLobatto[nGauss - 1] = 1.0;
            wGaussLobatto[nGauss] = 1.0;
        }
        else if (nGauss == 2)
        {
            xGaussLobatto[nGauss - 2] = -1.0;
            xGaussLobatto[nGauss - 1] = 0.0;
            xGaussLobatto[nGauss] = 1.0;

            wGaussLobatto[nGauss - 2] = 1.0/3.0;
            wGaussLobatto[nGauss - 1] = 4.0/3.0;
            wGaussLobatto[nGauss] = 1.0/3.0;
        }
    }

    void basisFc(double a, double b, int elemType)
    {
        switch (elemType)
        {
        case 3:
        {
            for (int i = 0; i <= mathVar::orderElem; i++)
            {
                if (i == 0)
                {
                    mathVar::B[i] = 1.0;
                }
                else if (i == 1)
                {
                    mathVar::B[i] = (3.0*b + 1.0) / 2.0;
                }
                else if (i == 2)
                {
                    mathVar::B[i] = a * (1 - b);
                }
                else if (i == 3)
                {
                    mathVar::B[i] = -0.5 + b + 5.0*pow(b, 2) / 2.0;
                }
                else if (i == 4)
                {
                    mathVar::B[i] = a*(1-b)*(1.5+2.5*b);
                }
                else if (i == 5)
                {
                    mathVar::B[i] = (1.5*a*a-0.5)*(1-b)*(1-b);
                }
                else if (i == 6)
                {
                    mathVar::B[i] = -(3.0/8.0)-(15.0/8.0)*b+(15.0/8.0)*b*b+(35.0/8.0)*b*b*b;
                }
                else if (i == 7)
                {
                    mathVar::B[i] = a*(1-b)*(0.25+4.5*b+(21.0/4.0)*b*b);
                }
                else if (i == 8)
                {
                    mathVar::B[i] = (1.5*a*a-0.5*(1-b)*(1-b))*(2.5+3.5*b);
                }
                else if (i == 9)
                {
                    mathVar::B[i] = (2.5*a*a*a-1.5)*pow((1-b),3.0);
                }
            }
        }
        break;
        case 4:
        {
            for (int i = 0; i <= mathVar::orderElem; i++)
            {
                if (i == 0)
                {
                    mathVar::B[i] = 1.0;
                }
                else if (i == 1)
                {
                    mathVar::B[i] = a;
                }
                else if (i == 2)
                {
                    mathVar::B[i] = b;
                }
                else if (i == 3)
                {
                    mathVar::B[i] = a * b;
                }
            }
        }
        break;
        default:
            break;
        }
    }

    void dBasisFc(double a, double b, int elemType)
    {
        switch (elemType)
        {
        case 3:
        {
            for (int i = 0; i <= mathVar::orderElem; i++)
            {
                if (i == 0)
                {
                    mathVar::dBa[i] = 0;
                    mathVar::dBb[i] = 0;
                }
                else if (i == 1)
                {
                    mathVar::dBa[i] = 0;
                    mathVar::dBb[i] = 3.0 / 2.0;
                }
                else if (i == 2)
                {
                    mathVar::dBa[i] = 1 - b;
                    mathVar::dBb[i] = -a;
                }
                else if (i == 3)
                {
                    mathVar::dBa[i] = 0;
                    mathVar::dBb[i] = 1+5*b;
                }
                else if (i == 4)
                {
                    mathVar::dBa[i] = 1.5+b-2.5*b*b;
                    mathVar::dBb[i] = a*(1-5*b);
                }
                else if (i == 5)
                {
                    mathVar::dBa[i] = 3*a*(1-b)*(1-b);
                    mathVar::dBb[i] = 2*(b-1)*(1.5*a*a-0.5);
                }
                else if (i == 6)
                {
                    mathVar::dBa[i] = 0;
                    mathVar::dBb[i] = (-15.0/8.0)+(15.0/4.0)*b+(35.0/24.0)*b*b;
                }
                else if (i == 7)
                {
                    mathVar::dBa[i] = (1-b)*(0.25+4.5*b+(21.0/4.0)*b*b);
                    mathVar::dBb[i] = a*(17.0/4.0+3*b/8.0+7*b*b/4.0);
                }
                else if (i == 8)
                {
                    mathVar::dBa[i] = 3*a*(2.5+3.5*b);
                    mathVar::dBb[i] = (21.0/4.0)*(a*a-b*b)+9.0*b/2.0+3.0/4.0;
                }
                else if (i == 9)
                {
                    mathVar::dBa[i] = (15.0/2.0)*a*a*pow(1-b,3.0);
                    mathVar::dBb[i] = -3*pow(1-b,2.0)*(5*a*a*a/2.0-1.5);
                }
            }
        }
        break;
        case 4:
        {
            for (int i = 0; i <= mathVar::orderElem; i++)
            {
                if (i == 0)
                {
                    mathVar::dBa[i] = 0.0;
                    mathVar::dBb[i] = 0.0;
                }
                else if (i == 1)
                {
                    mathVar::dBa[i] = 1.0;
                    mathVar::dBb[i] = 0.0;
                }
                else if (i == 2)
                {
                    mathVar::dBa[i] = 0.0;
                    mathVar::dBb[i] = 1.0;
                }
                else if (i == 3)
                {
                    mathVar::dBa[i] = b;
                    mathVar::dBb[i] = a;
                }
            }
        }
        break;
        default:
            break;
        }
    }

    double J2DCal(int elem, double a, double b)
    {
        double Jacobi(0.0), xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);
        int elemType(auxUlti::checkType(elem));
        if (elemType == 4)  //Quad
        {
            std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
            std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
            std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);
            std::tie(xD, yD) = auxUlti::getElemCornerCoord(elem, 3);

            Jacobi = math::jacobianQuad(xA, xB, xC, xD, yA, yB, yC, yD, a, b);
        }
        else if (elemType == 3)  //Tri
        {
            std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
            std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
            std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);

            Jacobi = math::jacobianTri(xA, xB, xC, yA, yB, yC, a, b);
        }
        return Jacobi;
    }

    std::tuple<double, double> J1DCal(int edge)
    {
        int master(0), servant(0), masterIndex(0), servantIndex(0), masterType(0), servantType(0);
        double JMaster(0.0), JServant(0.0), xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);

        if (meshVar::ineled[edge][0]>meshVar::ineled[edge][1])
        {
            master = meshVar::ineled[edge][0];
            servant = meshVar::ineled[edge][1];
        }
        else if (meshVar::ineled[edge][0]<meshVar::ineled[edge][1])
        {
            master = meshVar::ineled[edge][1];
            servant = meshVar::ineled[edge][0];
        }

        masterIndex = auxUlti::findEdgeOrder(master, edge);
        masterType = auxUlti::checkType(master);

        if (masterType==4)
        {
            std::tie(xA, yA) = auxUlti::getElemCornerCoord(master, 0);
            std::tie(xB, yB) = auxUlti::getElemCornerCoord(master, 1);
            std::tie(xC, yC) = auxUlti::getElemCornerCoord(master, 2);
            std::tie(xD, yD) = auxUlti::getElemCornerCoord(master, 3);

            JMaster = math::jacobian1DQuad(masterIndex, xA, xB, xC, xD, yA, yB, yC, yD);
        }
        else if (masterType==3)
        {
            std::tie(xA, yA) = auxUlti::getElemCornerCoord(master, 0);
            std::tie(xB, yB) = auxUlti::getElemCornerCoord(master, 1);
            std::tie(xC, yC) = auxUlti::getElemCornerCoord(master, 2);

            JMaster = math::jacobian1DTri(masterIndex, xA, xB, xC, yA, yB, yC);
        }

        if (servant >= 0)
        {
            servantIndex = auxUlti::findEdgeOrder(servant, edge);
            servantType = auxUlti::checkType(servant);
            if (servantType == 4)
            {
                std::tie(xA, yA) = auxUlti::getElemCornerCoord(servant, 0);
                std::tie(xB, yB) = auxUlti::getElemCornerCoord(servant, 1);
                std::tie(xC, yC) = auxUlti::getElemCornerCoord(servant, 2);
                std::tie(xD, yD) = auxUlti::getElemCornerCoord(servant, 3);

                JServant = math::jacobian1DQuad(servantIndex, xA, xB, xC, xD, yA, yB, yC, yD);
            }
            else if (servantType == 3)
            {
                std::tie(xA, yA) = auxUlti::getElemCornerCoord(servant, 0);
                std::tie(xB, yB) = auxUlti::getElemCornerCoord(servant, 1);
                std::tie(xC, yC) = auxUlti::getElemCornerCoord(servant, 2);

                JServant = math::jacobian1DTri(servantIndex, xA, xB, xC, yA, yB, yC);
            }
        }
        else
        {
            JServant = 0.0;
        }

        return std::make_tuple(JMaster, JServant);
    }

    double iniIntegral(int elemType, int order)
    {
        double w1(0.0), w2(0.0), integral(0.0);
        switch (elemType)
        {
        case 3:
        {
            for (int na = 0; na <= mathVar::nGauss2D; na++)
            {
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                {
                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                    w1 = mathVar::wGaussPts[nanb][0];
                    w2 = mathVar::wGaussPts[nanb][1];
                    integral += w1 * w2* mathVar::BPts_Tri[order][nanb];
                }
            }
        }
        break;
        case 4:
        {
            for (int na = 0; na <= mathVar::nGauss2D; na++)
            {
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                {
                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                    w1 = mathVar::wGaussPts[nanb][0];
                    w2 = mathVar::wGaussPts[nanb][1];
                    integral += w1 * w2* mathVar::BPts_Quad[order][nanb];
                }
            }
        }
        break;
        default:
            break;
        }
        return integral;
    }

    double volumeInte(std::vector< std::vector<double> > &Fvalue, int elem)
    {
        double J2D(0.0), w1(0.0), w2(0.0), integral(0.0);
        for (int na = 0; na <= mathVar::nGauss2D; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
            {
                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                w1 = mathVar::wGaussPts[nanb][0];
                w2 = mathVar::wGaussPts[nanb][1];
                J2D = meshVar::J2D[elem][nanb];
                integral += w1 * w2*(J2D)*Fvalue[na][nb];
            }
        }
        return integral;
    }

    double jacobianQuad(double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD, double a, double b)
    {
        double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
        double jQuad(0.0);

        dxa = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*b;
        dxb = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*a;
        dya = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*b;
        dyb = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*a;
        jQuad = dxa * dyb - dxb * dya;
        return jQuad;
    }

    double jacobianTri(double xA, double xB, double xC, double yA, double yB, double yC, double a, double b)
    {
        double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
        double jTri(0.0);

        dxa = (1 - b)*(xB - xA) / 4.0;
        dxb = a * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
        dya = (1 - b)*(yB - yA) / 4.0;
        dyb = a * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
        jTri = dxa * dyb - dxb * dya;
        return jTri;
    }

    double jacobian1DQuad(int edgeIndex, double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD)
    {
        double dx(0.0), dy(0.0), C(0.0);
        double jacobi(0.0);
        if ((edgeIndex == 0) || (edgeIndex == 2))  //AB or CD
        {
            if (edgeIndex == 0)
            {
                C = -1.0;
            }
            else if (edgeIndex == 2)
            {
                C = 1.0;
            }
            dx = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*C;
            dy = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*C;
        }
        else if ((edgeIndex == 1) || (edgeIndex == 3))  //BC or DA
        {
            if (edgeIndex == 1)
            {
                C = 1.0;
            }
            else if (edgeIndex == 3)
            {
                C = -1.0;
            }
            dx = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*C;
            dy = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*C;
        }

        jacobi = std::sqrt(pow(dx, 2) + pow(dy, 2));
        return jacobi;
    }

    double jacobian1DTri(int edgeIndex, double xA, double xB, double xC, double yA, double yB, double yC)
    {
        double dx(0.0), dy(0.0), C(0.0);
        double jacobi(0.0);
        if ((edgeIndex == 1) || (edgeIndex == 2))  //BC or CA
        {
            if (edgeIndex == 1)
            {
                C = 1.0;
            }
            else if (edgeIndex == 2)
            {
                C = -1.0;
            }
            dx = C * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
            dy = C * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
        }
        else if (edgeIndex == 0)  //AB
        {
            C = -1.0;
            dx = (1 - C)*(xB - xA) / 4.0;
            dy = (1 - C)*(yB - yA) / 4.0;
        }

        jacobi = std::sqrt(pow(dx, 2) + pow(dy, 2));
        return jacobi;
    }

    double calcMeanPriVar(int element, int varId)
    {
        /* Ham tinh gia tri trung binh cua bien primary tren cell (trung binh gia tri tai cac Gauss point)
         * Hien tai chi co bien T, varId = 1*/
        double sum(0.0);
        int count(0);
        switch (varId)
        {
        case 1:
        {
            for (int na = 0; na <= mathVar::nGauss2D; na++)
            {
                for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
                {
                    int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                    sum+=volumeFields::T[element][nanb];
                    count++;
                }
            }
            break;
        }
        case 2:
        {
            break;
        }
        default:
            break;
        }

        return (sum/count);
    }

    std::vector<double> SolveSysEqs(std::vector< std::vector<double> > &a, std::vector<double> &b)
    {
        int n = static_cast<int>(b.size());
        double eMax(1e-9), e(1.0), sum(0.0), xi(0.0);
        std::vector<double> results(n, 1.0);
        int counter(0);
        //eMax = fabs(*std::min_element(b.begin(), b.end())) / 1e3;

        while (e>eMax && counter<=15)
        {
            for (int i = 0; i < n; i++)
            {
                sum = b[i];
                for (int j = 0; j < n; j++)
                {
                    if (i != j)
                    {
                        sum -= a[i][j] * results[j];
                    }
                }
                xi = sum / a[i][i];
                results[i] = xi;
            }
            e = math::errorGS(a, b, results);
            counter++;
        }
        return results;
    }

    std::tuple<double,double> solveSys2Eqs(double *A, double *B)
    {
        /* Form |A[0] A[1]| |x1| = |B[0]|
         *      |A[2] A[3]| |x2| = |B[1]|
        */
        double x2=(B[0]*A[2]-A[0]*B[1])/(A[2]*A[1]-A[0]*A[3]);
        double x1=(B[0]-A[1]*x2)/A[0];
        return std::make_tuple(x1,x2);
    }

    double errorGS(std::vector< std::vector<double> > &a, std::vector<double> &b, std::vector<double> &No)
    {
        /* Calculate a*No-b */
        int n = static_cast<int>(b.size());
        double error(1.0);
        std::vector<double> R(n, 0.0);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                R[i] += a[i][j] * No[j];
            }
            R[i] = fabs(R[i] - b[i]);
        }
        error = *std::max_element(R.begin(), R.end());  //find max value of vector
        return error;
    }

    double CalcVisCoef(double T)
    {
        double mu(0.0);
        if (material::viscousityModel::sutherland)
        {
            mu = material::viscosityCoeff::Sutherland::As*pow(T, 1.5) / (T + material::viscosityCoeff::Sutherland::Ts);
        }
        else if (material::viscousityModel::power_VHS)
        {
            double Avogadro(6.02214076e23); //Note, molMass is g/mol
            double muRef(
                        15*pow((3.14*(material::viscosityCoeff::powerLaw_VHS::molMass/(Avogadro*1000))*
                                material::viscosityCoeff::powerLaw_VHS::kBoltzmann*
                                material::viscosityCoeff::powerLaw_VHS::TRef),
                               0.5)
                        /
                        (2*3.14*
                         pow(material::viscosityCoeff::powerLaw_VHS::dRef,2)*
                         (5-2*material::viscosityCoeff::powerLaw_VHS::omega)*
                         (7-2*material::viscosityCoeff::powerLaw_VHS::omega)
                         )
                        );

            mu = muRef*pow(T/material::viscosityCoeff::powerLaw_VHS::TRef,
                           material::viscosityCoeff::powerLaw_VHS::omega);
        }
        else if (material::viscousityModel::constant)
        {
            mu = material::viscosityCoeff::constant::mu;
        }

        if (mu!=mu)
        {
            std::cout << "Negative T = " << T << std::endl;
            exit(1);
        }

        return mu;
    }

    double CalcT(int elem, double a, double b)
    {
        double T(0.0), rhoGs(0.0), rhouGs(0.0), rhovGs(0.0), rhoEGs(0.0);
        rhoGs = math::pointValue(elem, a, b, 1, 2);
        rhouGs = math::pointValue(elem, a, b, 2, 2);
        rhovGs = math::pointValue(elem, a, b, 3, 2);
        rhoEGs = math::pointValue(elem, a, b, 4, 2);
        T = (material::gamma - 1)*(rhoEGs - 0.5*(pow(rhouGs, 2) + pow(rhovGs, 2)) / rhoGs) / (material::R*rhoGs);
        return T;
    }

    double CalcTFromConsvVar(double rho, double rhou, double rhov, double rhoE)
    {
        double out((material::gamma - 1)*(rhoE - 0.5*(pow(rhou, 2) + pow(rhov, 2)) / rho) / (material::R*rho));
        if (out < limitVal::TDwn)
        {
            //std::cout<<rho<<", "<<rhou<<", "<<rhov<<", "<<rhoE<< std::endl;
            //std::cout << "Negative T=" << out<< std::endl;
            //exit(1);
            out=limitVal::TDwn;
        }
        else if (out > limitVal::TUp)
        {
            out=limitVal::TUp;
        }
        return out;
    }

    double CalcP(double T, double rho)
    {
        double p = T * (material::R*rho);
        return p;
    }

    std::tuple<double, double, double, double> Calc_dxydab(int elem, double a, double b)
    {
        double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
        double xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);
        int elemType(0);

        elemType = auxUlti::checkType(elem);
        std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
        std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
        std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);

        if (elemType==4)  //Quad
        {
            std::tie(xD, yD) = auxUlti::getElemCornerCoord(elem, 3);
            dxa = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*b;
            dxb = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*a;
            dya = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*b;
            dyb = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*a;
        }
        else if (elemType == 3)  //Tri
        {
            dxa = (1 - b)*(xB - xA) / 4.0;
            dxb = a * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
            dya = (1 - b)*(yB - yA) / 4.0;
            dyb = a * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
        }

        return std::make_tuple(dxa, dxb, dya, dyb);
    }

    double Calc_dBxdBy(int elem, int order, int na, int nb, int opt)
    {
        int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
        double dB(0.0);
        int elemType(auxUlti::checkType(elem));
        switch (elemType)
        {
        case 3:
        {
            if (opt == 1)  //x direction
            {
                dB = (1 / meshVar::J2D[elem][nanb]) * (mathVar::dBaPts_Tri[order][nanb] * meshVar::dyb[elem][nanb] - mathVar::dBbPts_Tri[order][nanb] * meshVar::dya[elem][nanb]);
            }
            else if (opt == 2)  //y direction
            {
                dB = (1 / meshVar::J2D[elem][nanb]) * (mathVar::dBbPts_Tri[order][nanb] * meshVar::dxa[elem][nanb] - mathVar::dBaPts_Tri[order][nanb] * meshVar::dxb[elem][nanb]);
            }
        }
        break;
        case 4:
        {
            if (opt == 1)  //x direction
            {
                dB = (1 / meshVar::J2D[elem][nanb]) * (mathVar::dBaPts_Quad[order][nanb] * meshVar::dyb[elem][nanb] - mathVar::dBbPts_Quad[order][nanb] * meshVar::dya[elem][nanb]);
            }
            else if (opt == 2)  //y direction
            {
                dB = (1 / meshVar::J2D[elem][nanb]) * (mathVar::dBbPts_Quad[order][nanb] * meshVar::dxa[elem][nanb] - mathVar::dBaPts_Quad[order][nanb] * meshVar::dxb[elem][nanb]);
            }
        }
        break;
        default:
            break;
        }
        return dB;
    }

    double surfaceInte(std::vector<double> &Fvalue, int edge)
    {
        double inte(0.0), J(meshVar::J1D[edge]), w(0.0);
        for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
        {
            w = mathVar::wGaussSur[nG];
            inte += w * Fvalue[nG] * J;
        }
        return inte;
    }

    double pointValue(int element, double a, double b, int valType, int valKind)
    {
        //Compute primary variables from conservative variables
        double out(0.0);
        if (valKind==1)  //primary variables
        {
            if (valType == 1)  //rho
            {
                out= math::calcConsvVarWthLimiter(element, a, b, valType);
            }
            else if (valType == 2)  //u
            {
                double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
                    rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2));
                out = rhouVal / rhoVal;
            }
            else if (valType == 3)  //v
            {
                double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
                    rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3));
                out = rhovVal / rhoVal;
            }
            else if (valType == 4)  //e
            {
                double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
                    rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
                    rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
                    rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));
                    out = material::Cv*math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
            }
            else if (valType == 5)  //p
            {
                double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
                    rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
                    rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
                    rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));
                    out = math::CalcP(CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal), rhoVal);
            }
            else if (valType == 6)  //T
            {
                double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
                    rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
                    rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
                    rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));

                    out = CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
            }
            else if (valType == 7)  //mu
            {
                double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
                    rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
                    rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
                    rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));
                double TVal(CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));

                out = math::CalcVisCoef(TVal);
                if (out < 0 || out != out)
                {
                    std::cout << "unphysical mu at cell " << element + meshVar::nelem1D + 1 << std::endl;
                    exit(1);
                }
            }
        }
        else if (valKind==2)  //conservative variables
        {
            out = math::calcConsvVarWthLimiter(element, a, b, valType);
        }
        return out;
    }

    double pointValueNoLimiter(int element, double a, double b, int valType)
    {
        double out(0.0);
        std::vector<double> Value(mathVar::orderElem + 1, 0.0);
        Value = auxUlti::getElementConserValuesOfOrder(element, valType);

        math::basisFc(a, b, auxUlti::checkType(element));
        for (int order = 0; order <= mathVar::orderElem; order++)
        {
            out += Value[order] * mathVar::B[order];
        }

        return out;
    }

    double vectorDotProduct(const std::vector<double> &a, const std::vector<double> &b)
    {
        double out(0.0);
        out = a[0] * b[0] + a[1] * b[1];
        return out;
    }

    std::tuple <double, double> internalSurfaceValue(int edge, int element, int nG, int valType, int valKind)
    {
        int masterElem(0), servantElem(0);
        std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
        double valPlus(0.0), valMinus(0.0), aMaster(0.0), bMaster(0.0), aServant(0.0), bServant(0.0);

        std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
        std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, servantElem, nG);

        if (masterElem==element)  //considering element is master
        {
            valPlus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
            valMinus = math::pointValue(servantElem, aServant, bServant, valType, valKind);
        }
        else
        {
            valMinus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
            valPlus = math::pointValue(servantElem, aServant, bServant, valType, valKind);
        }

        return std::make_tuple(valPlus, valMinus);
    }

    std::tuple <double, double> internalSurfaceDerivativeValue(int edge, int element, int nG, int valType, int dir)
    {
        /* ValType la id cua bien can tinh:
         * 1: (mu)dRho
         * 2: (mu)dRhou
         * 3: (mu)dRhov
         * 4: (mu)dRhoE
         * 5: dRho --> day la bien phu khi mass diffusion on: Sm = div(rho)
        */

        int masterElem(0), servantElem(0);
        std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
        double valPlus(0.0), valMinus(0.0), aMaster(0.0), bMaster(0.0), aServant(0.0), bServant(0.0);

        std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
        std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, servantElem, nG);

        if (masterElem == element)  //considering element is master
        {
            if (systemVar::auxVariables==1)
            {
                valPlus = math::pointAuxValue(masterElem, aMaster, bMaster, valType, dir);
                valMinus = math::pointAuxValue(servantElem, aServant, bServant, valType, dir);
            }
            else if (systemVar::auxVariables==2)
            {
                valPlus = math::BR2Fncs::pointAuxValue_sur(edge,masterElem,aMaster,bMaster,valType,dir);
                valMinus = math::BR2Fncs::pointAuxValue_sur(edge,servantElem, aServant, bServant, valType, dir);
            }
        }
        else
        {
            if (systemVar::auxVariables==1)
            {
                valMinus = math::pointAuxValue(masterElem, aMaster, bMaster, valType, dir);
                valPlus = math::pointAuxValue(servantElem, aServant, bServant, valType, dir);
            }
            else if (systemVar::auxVariables==2)
            {
                valMinus = math::BR2Fncs::pointAuxValue_sur(edge,masterElem,aMaster,bMaster,valType,dir);
                valPlus = math::BR2Fncs::pointAuxValue_sur(edge,servantElem, aServant, bServant, valType, dir);
            }
        }

        return std::make_tuple(valPlus, valMinus);
    }

    double SurfaceValueFromMaster(int edge, int nG, int valType, int valKind)
    {
        int masterElem(0), servantElem(0);
        std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
        double valPlus(0.0), aMaster(0.0), bMaster(0.0);
        std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
        valPlus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
        return valPlus;
    }

    double CalcSpeedOfSound(double T)
    {
        double C(sqrt(material::gamma*material::R*T));
        return C;
    }

    std::vector<double> vectorSum(std::vector<double> &a, std::vector<double> &b)
    {
        static int size(static_cast<int>(a.size()));
        std::vector<double> out(size, 0.0);
        for (int i = 0; i < size; i++)
        {
            out[i] = a[i] + b[i];
        }
        return out;
    }

    double calcRhouvEDeriv(double dRhoUVal, double dRhoVal, double rhoUVal, double rhoVal)
    {
        /*NOTE! direction of derivatives based on input values, for example:
        d(u)/dx = [d(rhou)/dx - d(rho)/dx.(rhou)/rho]/(rho)
        so input values are
        d(rho)/dx	----> dRhoVal
        d(rhou)/dx	----> dRhoUVal
        rhou		----> rhoUVal
        rho			----> rhoVal*/
        double outVal(0.0);
        outVal = (dRhoUVal - dRhoVal * rhoUVal / rhoVal) / rhoVal;
        return outVal;
    }

    double calcTDeriv(double dE, double du, double dv, double u, double v)
    {
        /*NOTE! direction of derivatives based on input values, for example:
        d(T)/dx = f(dE/dx, du/dx, dv/dx, u, v)
        so input values are
        dE/dx	----> dE
        du/dx	----> du
        dv/dx	----> dv*/
        double outVal(0.0);
        outVal = (material::gamma - 1)*(dE - (u*du + v * dv)) / material::R;
        return outVal;
    }

    /**
     * @brief Function calculates auxiliary variables at abitrary point.
     * @param element: element Id
     * @param a: Point's coordinates
     * @param b: Point's coordinates
     * @param valType: id of variables
     * @param dir: direction, 1 for Ox, 2 for Oy
     *
     * Id	|keyWord	|
     * -----|-----------|
     * 1    |(mu)dRho	|
     * 2    |(mu)dRhou	|
     * 3    |(mu)dRhov	|
     * 4    |(mu)dRhoE	|
     * 5    |dRho		|
     *
     * @return
     */
    double pointAuxValue(int element, double a, double b, int valType, int dir)
    {
        /* ValType la id cua bien can tinh:
         * 1: (mu)dRho
         * 2: (mu)dRhou
         * 3: (mu)dRhov
         * 4: (mu)dRhoE
        */

        double out(0.0);
        std::vector<double> Value(mathVar::orderElem + 1, 0.0);

        Value = auxUlti::getElementAuxValuesOfOrder(element, valType, dir);

        math::basisFc(a, b, auxUlti::checkType(element));

        //Limiter
        double theta(1.0);
        if (flowProperties::massDiffusion)
        {
            theta=theta1Arr[element]*theta2Arr[element];
        }

        for (int order = 1; order <= mathVar::orderElem; order++)
        {
            out += Value[order] * mathVar::B[order]*theta;
        }
        out+=Value[0];

        return out;
    }

    double pointDerivRho(int element, double a, double b, int dir)
    {
        double out(0.0);
        std::vector<double> Value(mathVar::orderElem + 1, 0.0);

        if (dir==1)  //Ox direction
        {
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
            {
                Value[iorder] = BR1Vars::rhoX[element][iorder];
            }
        }
        else if (dir==2)  //Ox direction
        {
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
            {
                Value[iorder] = BR1Vars::rhoY[element][iorder];
            }
        }

        math::basisFc(a, b, auxUlti::checkType(element));
        //Limiter
        double theta(1.0);
        if (flowProperties::massDiffusion)
        {
            theta=theta1Arr[element]*theta2Arr[element];
        }

        for (int order = 1; order <= mathVar::orderElem; order++)
        {
            out += Value[order] * mathVar::B[order]*theta;
        }
        out += Value[0];

        return out;
    }

    /**
     * @brief Function calculates thermal conductivity
     * @param muVal: dynamic viscosity \f$\mu\f$
     * @return
     */
    double calcThermalConductivity(double muVal)
    {
        double k(0.0);
        k = material::Cp*muVal / material::Pr;
        return k;
    }

    /**
     * @brief Function calculates devivative of T from derivative of conservative variables.
     * @param dRhoE: div of total energy \f$\nabla (\rho E)\f$
     * @param dRho: div of density \f$\nabla (\rho)\f$
     * @param rhoE: total energy \f$\rho E\f$
     * @param rho: density \f$\rho\f$
     * @return
     */
    double calcDivTFromAuxVariable(double dRhoE, double dRho, double rhoE, double rho)
    {
        //Dung cach cua Thay, coi u^2 la constant. Vi dang tinh cho bien phu la S = mu*div(U) nen gia tri nay ma mu*divT
        return ((dRhoE-dRho*rhoE/rho)/(rho*material::Cv));
    }

    /**
     * @brief Function calculate Mean Free Path
     * @param mu: dynamic viscosity \f$\mu\f$
     * @param rho: density \f$\rho\f$
     * @param T: temperature T
     * @return
     */
    double calcMeanFreePath(double mu, double rho, double T)
    {
        return (mu*pow(3.1416/(2*material::R*T),0.5)/rho);
    }

    /**
     * @brief Function solve quadratic equation.
     *
     * Form of equation: \f$Ax^2 + Bx + C = 0\f$
     *
     * @param A: coefficient
     * @param B: coefficient
     * @param C: coefficient
     * @return
     */
    std::tuple<bool, double, double> solvQuadraticEq(double A, double B, double C)
    {
        double delta(B*B - 4 * A*C), root1(0.0), root2(0.0);
        bool realRoot(true);

        if (A != 0.0)
        {
            if (delta>0)
            {
                root1 = ((-B + sqrt(delta)) / (2 * A));
                root2 = ((-B - sqrt(delta)) / (2 * A));
            }
            else if (delta == 0.0)
            {
                root1 = (-B / (2 * A));
                root2 = root1;
            }
            else
            {
                realRoot = false;
            }
        }
        else if (A == 0.0 && B !=0.0)
        {
            root1 = -C / B;
            root2 = root1;
        }
        else
        {
            realRoot = false;
        }
        return std::make_tuple(realRoot, root1, root2);
    }

    double centerValue(int element, int valType, int valKind)
    {
        double xC(-1.0 / 3.0), yC(1.0 / 3.0), output(0.0);

        if (auxUlti::checkType(element) == 4)
        {
            xC = 0.0;
            yC = 0.0;
        }
        output = math::pointValue(element, xC, yC, valType, valKind);
        return output;
    }

    /**
     * @brief Function calculates value of auxiliary variable at cell's center.
     * @param element: element Id
     * @param valType: id of variables
     * @param dir: direction, 1 for Ox, 2 for Oy
     * Id	|keyWord	|
     * -----|-----------|
     * 1    |(mu)dRho	|
     * 2    |(mu)dRhou	|
     * 3    |(mu)dRhov	|
     * 4    |(mu)dRhoE	|
     * 5    |dRho		|
     * @return
     */
    double centerAuxValue(int element, int valType, int dir)
    {
        double xC(-1.0 / 3.0), yC(1.0 / 3.0), output(0.0);

        if (auxUlti::checkType(element) == 4)
        {
            xC = 0.0;
            yC = 0.0;
        }
        if (systemVar::auxVariables==1)
        {
            output = math::pointAuxValue(element, xC, yC, valType, dir);
        }
        else if (systemVar::auxVariables==2)
        {
            output = math::BR2Fncs::pointAuxValue_vol(element, xC, yC, valType, dir);
        }
        return output;
    }

    /**
     * @brief Function maps a point in standard coordinates system to real coordinates system.
     * @param element: element Id
     * @param aCoor: point's coordinate in standard coordinates system
     * @param bCoor: point's coordinate in standard coordinates system
     * @return
     */
    std::tuple<double, double> directMapping(int element, double aCoor, double bCoor)
    {
        int elemType(auxUlti::checkType(element));
        double xCoor(0.0), yCoor(0.0), xA(0.0), xB(0.0), xC(0.0),
            yA(0.0), yB(0.0), yC(0.0);

        std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
        std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
        std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

        if (elemType == 3) //Tri element
        {
            xCoor = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.5*(1 + bCoor)*xC;
            yCoor = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.5*(1 + bCoor)*yC;
        }
        else if (elemType == 4) //Quad element
        {
            double xD(0.0), yD(0.0);
            std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

            xCoor = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.25*(1 - aCoor)*(1 + bCoor)*xD + 0.25*(1 + aCoor)*(1 + bCoor)*xC;
            yCoor = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.25*(1 - aCoor)*(1 + bCoor)*yD + 0.25*(1 + aCoor)*(1 + bCoor)*yC;
        }
        return std::make_tuple(xCoor, yCoor);
    }

    /**
     * @brief Function maps a point in real coordinates system to standard coordinates system (old version).
     * @param edge: edge Id
     * @param element: element Id
     * @param xCoor: point's coordinate in real coordinates system
     * @param yCoor: point's coordinate in real coordinates system
     * @return
     */
    std::tuple<double, double> mappingRealToStd(int edge, int element, double xCoor, double yCoor)
    {
        double xA(0.0), xB(0.0), xC(0.0),
            yA(0.0), yB(0.0), yC(0.0),
            aCoor(0.0), bCoor(0.0),
            xCal(0.0), yCal(0.0), eX(100.0), eY(100.0);
        int elemType(auxUlti::checkType(element));
        std::vector<std::vector<double>> vectorGaussPoints(mathVar::nGauss1D + 1, std::vector<double>(2, 0.0));

        std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
        std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
        std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

        vectorGaussPoints = auxUlti::getVectorGaussSurfCoor(edge, element);

        if (elemType==3) //tri
        {
            for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
            {
                aCoor = vectorGaussPoints[nG][0];
                bCoor = vectorGaussPoints[nG][1];
                xCal = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.5*(1 + bCoor)*xC;
                yCal = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.5*(1 + bCoor)*yC;
                eX = fabs(xCal - xCoor) * 100 / xCoor;
                eY = fabs(yCal - yCoor) * 100 / yCoor;
                if ((eX < 0.005) && (eY < 0.005))
                {
                    break;
                }
            }
        }
        else if (elemType==4) //quad
        {
            double xD(0.0), yD(0.0);
            std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

            for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
            {
                aCoor = vectorGaussPoints[nG][0];
                bCoor = vectorGaussPoints[nG][1];
                xCal = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.25*(1 - aCoor)*(1 + bCoor)*xD + 0.25*(1 + aCoor)*(1 + bCoor)*xC;
                yCal = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.25*(1 - aCoor)*(1 + bCoor)*yD + 0.25*(1 + aCoor)*(1 + bCoor)*yC;
                eX = fabs((xCal - xCoor) * 100 / xCoor);
                eY = fabs((yCal - yCoor) * 100 / yCoor);
                if ((eX < 0.05) && (eY < 0.05))
                {
                    break;
                }
            }
        }

        return std::make_tuple(aCoor, bCoor);
    }

    /**
     * @brief Function maps a point in real coordinates system to standard coordinates system (current version).
     * @param element: element Id
     * @param xCoor: point's coordinate in real coordinates system
     * @param yCoor: point's coordinate in real coordinates system
     * @return
     */
    std::tuple<double, double> inverseMapping(int element, double xCoor, double yCoor)
    {
        int elemType(auxUlti::checkType(element));
        double aCoor(0.0), bCoor(0.0), xA(0.0), xB(0.0), xC(0.0),
            yA(0.0), yB(0.0), yC(0.0),
            Ax(0.0), Bx(0.0), Cx(0.0), Dx(0.0),
            Ay(0.0), By(0.0), Cy(0.0), Dy(0.0), AA(0.0), BB(0.0), CC(0.0);
        bool realRoot(true);

        std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
        std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
        std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

        if (elemType == 3) //Tri element
        {
            Ax = (xA - xB) / 4.0;
            Bx = (xB - xA) / 4.0;
            Cx = (2 * xC - xA - xB) / 4.0;
            Dx = (2 * xC + xA + xB) / 4.0 - xCoor;
            Ay = (yA - yB) / 4.0;
            By = (yB - yA) / 4.0;
            Cy = (2 * yC - yA - yB) / 4.0;
            Dy = (2 * yC + yA + yB) / 4.0 - yCoor;
        }
        else if (elemType == 4) //Quad element
        {
            double xD(0.0), yD(0.0);
            std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

            Ax = (xA - xB - xD + xC) / 4.0;
            Bx = (- xA + xB - xD + xC) / 4.0;
            Cx = (- xA - xB + xD + xC) / 4.0;
            Dx = (xA + xB + xD + xC) / 4.0 - xCoor;
            Ay = (yA - yB - yD + yC) / 4.0;
            By = (-yA + yB - yD + yC) / 4.0;
            Cy = (-yA - yB + yD + yC) / 4.0;
            Dy = (yA + yB + yD + yC) / 4.0 - yCoor;
        }

        //solve a
        AA = -Ax * By + Ay * Bx;
        BB = -Ax * Dy + Bx * Cy - Cx * By + Ay * Dx;
        CC = Dx * Cy - Cx * Dy;
        if (fabs(AA) < 1e-12 && fabs(BB) < 1e-12 && fabs(CC) < 1e-12)
        {
            aCoor = -1;
            bCoor = 1;
        }
        else
        {
            std::vector<double> vectorRoot_a(2, 0.0);
            std::tie(realRoot, vectorRoot_a[0], vectorRoot_a[1]) = math::solvQuadraticEq(AA, BB, CC);
            if (realRoot)
            {
                for (int i = 0; i < 2; i++)
                {
                    double error(fabs(fabs(vectorRoot_a[i]) - 1) * 100);
                    if ((fabs(vectorRoot_a[i]) <= 1) || (error <= 0.00001))
                    {
                        aCoor = vectorRoot_a[i];
                        break;
                    }
                }
            }
            else
            {
                /*
                std::cout << "Element " << element << std::endl;
                std::cout << "A: " << AA << " B: " << BB << " C: " << CC << std::endl;
                std::cout << "x: " << xCoor << " y: " << yCoor << std::endl;
                */
                std::string str("mapping error occured_a");
                exitDG(str);
            }

            if (fabs(aCoor*Ay + Cy) <= 1e-12)
            {
                //std::cout << "aCoor: " << aCoor << " criteria: " << fabs((aCoor*By + Dy) - (aCoor*Ay + Cy)) << std::endl;
                //bCoor = 2;
                AA = -Ay * Bx + Ax * By;
                BB = -Ay * Dx + By * Cx - Cy * Bx + Ax * Dy;
                CC = Dy * Cx - Cy * Dx;
                if (fabs(AA) < 1e-12 && fabs(BB) < 1e-12 && fabs(CC) < 1e-12)
                {
                    aCoor = -1;
                    bCoor = 1;
                }
                else
                {
                    std::vector<double> vectorRoot_a(2, 0.0);
                    std::tie(realRoot, vectorRoot_a[0], vectorRoot_a[1]) = math::solvQuadraticEq(AA, BB, CC);
                    if (realRoot)
                    {
                        for (int i = 0; i < 2; i++)
                        {
                            double error(fabs(fabs(vectorRoot_a[i]) - 1) * 100);
                            if ((fabs(vectorRoot_a[i]) <= 1) || (error <= 0.00001))
                            {
                                aCoor = vectorRoot_a[i];
                                break;
                            }
                        }
                    }
                    else
                    {
                        /*
                        std::cout << "Element " << element << std::endl;
                        std::cout << "A: " << AA << " B: " << BB << " C: " << CC << std::endl;
                        std::cout << "x: " << xCoor << " y: " << yCoor << std::endl;
                        */
                        std::string str("mapping error occured");
                        exitDG(str);
                    }

                    if (fabs(aCoor*Ax + Cx) == 0)
                    {
                        std::string str("mapping error occured");
                        exitDG(str);
                    }
                    else
                    {
                        bCoor = -(aCoor*Bx + Dx) / (aCoor*Ax + Cx);
                    }
                }
            }
            else
            {
                bCoor = -(aCoor*By + Dy) / (aCoor*Ay + Cy);
            }
        }

        return std::make_tuple(aCoor, bCoor);
    }

    std::tuple<double, double> inverseMapping_ForParallel(int edge, double xCoor, double yCoor)
    {
        int elemType(3), loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        double aCoor(0.0), bCoor(0.0), xA(0.0), xB(0.0), xC(0.0),
            yA(0.0), yB(0.0), yC(0.0),
            Ax(0.0), Bx(0.0), Cx(0.0), Dx(0.0),
            Ay(0.0), By(0.0), Cy(0.0), Dy(0.0), AA(0.0), BB(0.0), CC(0.0);
        bool realRoot(true);

        xA=parallelBuffer::xCoor[loc][0];
        xB=parallelBuffer::xCoor[loc][1];
        xC=parallelBuffer::xCoor[loc][2];
        yA=parallelBuffer::yCoor[loc][0];
        yB=parallelBuffer::yCoor[loc][1];
        yC=parallelBuffer::yCoor[loc][2];

        if (elemType == 3) //Tri element
        {
            Ax = (xA - xB) / 4.0;
            Bx = (xB - xA) / 4.0;
            Cx = (2 * xC - xA - xB) / 4.0;
            Dx = (2 * xC + xA + xB) / 4.0 - xCoor;
            Ay = (yA - yB) / 4.0;
            By = (yB - yA) / 4.0;
            Cy = (2 * yC - yA - yB) / 4.0;
            Dy = (2 * yC + yA + yB) / 4.0 - yCoor;
        }
        else if (elemType == 4) //Quad element
        {
            //PARALLEL RUNNING CHUA KIEM TRA DUOC CHO QUADRILATERAL ELEMENT
            double xD(0.0), yD(0.0);
            //std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

            Ax = (xA - xB - xD + xC) / 4.0;
            Bx = (- xA + xB - xD + xC) / 4.0;
            Cx = (- xA - xB + xD + xC) / 4.0;
            Dx = (xA + xB + xD + xC) / 4.0 - xCoor;
            Ay = (yA - yB - yD + yC) / 4.0;
            By = (-yA + yB - yD + yC) / 4.0;
            Cy = (-yA - yB + yD + yC) / 4.0;
            Dy = (yA + yB + yD + yC) / 4.0 - yCoor;
        }

        //solve a
        AA = -Ax * By + Ay * Bx;
        BB = -Ax * Dy + Bx * Cy - Cx * By + Ay * Dx;
        CC = Dx * Cy - Cx * Dy;
        if (fabs(AA) < 1e-12 && fabs(BB) < 1e-12 && fabs(CC) < 1e-12)
        {
            aCoor = -1;
            bCoor = 1;
        }
        else
        {
            std::vector<double> vectorRoot_a(2, 0.0);
            std::tie(realRoot, vectorRoot_a[0], vectorRoot_a[1]) = math::solvQuadraticEq(AA, BB, CC);
            if (realRoot)
            {
                for (int i = 0; i < 2; i++)
                {
                    double error(fabs(fabs(vectorRoot_a[i]) - 1) * 100);
                    if ((fabs(vectorRoot_a[i]) <= 1) || (error <= 0.00001))
                    {
                        aCoor = vectorRoot_a[i];
                        break;
                    }
                }
            }
            else
            {
                /*
                std::cout << "Element " << element << std::endl;
                std::cout << "A: " << AA << " B: " << BB << " C: " << CC << std::endl;
                std::cout << "x: " << xCoor << " y: " << yCoor << std::endl;
                */
                std::string str("mapping error occured_a");
                exitDG(str);
            }

            if (fabs(aCoor*Ay + Cy) <= 1e-12)
            {
                //std::cout << "aCoor: " << aCoor << " criteria: " << fabs((aCoor*By + Dy) - (aCoor*Ay + Cy)) << std::endl;
                //bCoor = 2;
                AA = -Ay * Bx + Ax * By;
                BB = -Ay * Dx + By * Cx - Cy * Bx + Ax * Dy;
                CC = Dy * Cx - Cy * Dx;
                if (fabs(AA) < 1e-12 && fabs(BB) < 1e-12 && fabs(CC) < 1e-12)
                {
                    aCoor = -1;
                    bCoor = 1;
                }
                else
                {
                    std::vector<double> vectorRoot_a(2, 0.0);
                    std::tie(realRoot, vectorRoot_a[0], vectorRoot_a[1]) = math::solvQuadraticEq(AA, BB, CC);
                    if (realRoot)
                    {
                        for (int i = 0; i < 2; i++)
                        {
                            double error(fabs(fabs(vectorRoot_a[i]) - 1) * 100);
                            if ((fabs(vectorRoot_a[i]) <= 1) || (error <= 0.00001))
                            {
                                aCoor = vectorRoot_a[i];
                                break;
                            }
                        }
                    }
                    else
                    {
                        /*
                        std::cout << "Element " << element << std::endl;
                        std::cout << "A: " << AA << " B: " << BB << " C: " << CC << std::endl;
                        std::cout << "x: " << xCoor << " y: " << yCoor << std::endl;
                        */
                        std::string str("mapping error occured");
                        exitDG(str);
                    }

                    if (fabs(aCoor*Ax + Cx) == 0)
                    {
                        std::string str("mapping error occured");
                        exitDG(str);
                    }
                    else
                    {
                        bCoor = -(aCoor*Bx + Dx) / (aCoor*Ax + Cx);
                    }
                }
            }
            else
            {
                bCoor = -(aCoor*By + Dy) / (aCoor*Ay + Cy);
            }
        }
        return std::make_tuple(aCoor, bCoor);
    }

    double vectorNorm(std::vector<double> vector)
    {
        int size(static_cast<int>(vector.size()));
        double norm(0.0);
        for (int i = 0; i < size; i++)
        {
            norm += pow(vector[i], 2);
        }
        norm = pow(norm, 0.5);
        return norm;
    }

    std::vector<double> traceOfMatrix(std::vector<std::vector<double>> &M)
    {
        int size(static_cast<int>(M.size()));
        std::vector<double> trace(size,0.0);
        for (int i=0; i<size;i++)
        {
            trace[i]=M[i][i];
        }
        return trace;
    }

    std::vector<std::vector<double>> transformationTensor(double nx, double ny)
    {
        /*
         * Transformation tensor S = I - nn
         * I: identity tensor
         * [ 1 0 0]
         * [ 0 1 0]
         * [ 0 0 1]
         * nn la outer product cua n va n = n*nT
         * nn = [nx] * [nx ny nz] = [nx*nx nx*ny nx*nz]
         *      [ny]                [ny*nx ny*ny ny*nz]
         *      [nz]                [nz*nx nz*ny nz*nz]
        */
        std::vector<std::vector<double>> S{{ 1-nx*nx, -nx*ny },
                                           { -nx*ny, 1-ny*ny }};
        return S;
    }

    std::vector<std::vector<double>> gradU(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
    {
        //Function calculates grad of U (is a vector field) from conservative vars and auxilary vars
        double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0);
        double drhox(0.0), drhoy(0.0),
            drhoux(0.0), drhouy(0.0),
            drhovx(0.0), drhovy(0.0);

        rhoVal = U[0];
        rhouVal = U[1];
        rhovVal = U[2];

        drhox = dUx[0];
        drhoy = dUy[0];

        drhoux = dUx[1];
        drhouy = dUy[1];

        drhovx = dUx[2];
        drhovy = dUy[2];

        dux = math::calcRhouvEDeriv(drhoux, drhox, rhouVal, rhoVal);
        duy = math::calcRhouvEDeriv(drhouy, drhoy, rhouVal, rhoVal);

        dvx = math::calcRhouvEDeriv(drhovx, drhox, rhovVal, rhoVal);
        dvy = math::calcRhouvEDeriv(drhovy, drhoy, rhovVal, rhoVal);

        std::vector<std::vector<double>> gradU{{dux, duy}, {dvx, dvy}};
        return gradU;
    }

    std::vector<std::vector<double>> transposeDoubleMatrix(std::vector<std::vector<double>> &M)
    {
        int rowNum(static_cast<int>(M.size())), colNum(static_cast<int>(M[0].size()));
        std::vector<std::vector<double>> MT(colNum,std::vector<double>(rowNum));
        for (int row=0; row<rowNum; row++)
        {
            for (int col=0; col<colNum; col++)
            {
                MT[col][row]=M[row][col];
            }
        }
        return MT;
    }

    namespace tensorMath {
        std::vector<double> tensorVectorDotProduct(std::vector<std::vector<double>>&S, std::vector<double>&V)
        {
            int size(static_cast<int>(V.size()));
            std::vector<double> Product(size,0.0);
            for (int i=0; i<size; i++)
            {
                for (int j=0; j<size; j++)
                {
                    Product[i]+=S[i][j]*V[j];
                }
            }
            return Product;
        }

        std::vector<double> vectorTensorDotProduct(std::vector<double>&V, std::vector<std::vector<double>>&S)
        {
            int size(static_cast<int>(V.size()));
            std::vector<double> Product(size,0.0);
            for (int i=0; i<size; i++)
            {
                for (int j=0; j<size; j++)
                {
                    Product[i]+=V[j]*S[j][i];
                }
            }
            return Product;
        }

        std::vector<std::vector<double>> sumOf2Tensor(std::vector<std::vector<double>>&S1, std::vector<std::vector<double>>&S2, int coefOfS1, int coefOfS2)
        {
            int rowNum(static_cast<int>(S1.size())), colNum(static_cast<int>(S1[0].size()));
            std::vector<std::vector<double>> sum(rowNum,std::vector<double>(colNum));
            for (int row=0; row<rowNum; row++)
            {
                for (int col=0; col<colNum; col++)
                {
                    sum[row][col]=coefOfS1*S1[row][col]+coefOfS2*S2[row][col];
                }
            }
            return sum;
        }

        std::vector<double> sumOf2Vector(std::vector<double>&V1, std::vector<double>&V2, int coefOfV1, int coefOfV2)
        {
            int rowNum(static_cast<int>(V1.size()));
            std::vector<double> sum(rowNum,0.0);
            for (int row=0; row<rowNum; row++)
            {
                sum[row]=coefOfV1*V1[row]+coefOfV2*V2[row];
            }
            return sum;
        }

        std::vector<std::vector<double>> numberTensorProduct(double num, std::vector<std::vector<double>>&S)
        {
            int rowNum(static_cast<int>(S.size())), colNum(static_cast<int>(S[0].size()));
            std::vector<std::vector<double>> product(rowNum,std::vector<double>(colNum));
            for (int row=0; row<rowNum; row++)
            {
                for (int col=0; col<colNum; col++)
                {
                    product[row][col]=S[row][col]*num;
                }
            }
            return product;
        }

        std::vector<double> numberVectorProduct(double num, std::vector<double>&V)
        {
            int rowNum(static_cast<int>(V.size()));
            std::vector<double> product(rowNum,0.0);
            for (int row=0; row<rowNum; row++)
            {
                product[row]=V[row]*num;
            }
            return product;
        }

        std::vector<std::vector<double>> unitTensor(int rowNum)
        {
            std::vector<std::vector<double>> I(rowNum,std::vector<double>(rowNum, 0.0));
            for (int row=0; row<rowNum; row++)
            {
                I[row][row]=1.0;
            }
            return I;
        }
    }

    namespace numericalFluxes
    {
        void calConvectiveFluxes(int edgeId, std::vector<double>&flux, std::vector<double> &fxP, std::vector<double> &fxM, std::vector<double> &fyP, std::vector<double> &fyM, std::vector<double> &UP, std::vector<double> &UM, double TP, double TM, std::vector<double> &vectorn)
        {
            /* Ham kiem soat buoc tinh convective flux
             * - P la phia ben trong, M la phia ben ngoai phan tu dang xet.
             * - Neu flux type la LxF va central, cac bien su dung duoc de o he toa do global:
             *      + flux          : output flux o he toa do global
             *      + fxP, fxM      : inviscid flux theo phuong Ox global
             *      + fyP, fyM      : inviscid flux theo phuong Oy global
             *      + UP, UM        : vector bien conservative o he toa do global
             *      + TP, TM        : nhiet do
             *      + vectorn       : vector phap tuyen theo he toa do global
             *
             * - Neu flux type la Roe hoac HLL, cac bien su dung duoc de o he toa do local:
             *      + flux          : output flux o he toa do global
             *      + fxP, fxM      : inviscid flux theo phuong Ox local (phuong normal)
             *      + fyP, fyM      : --------
             *      + UP, UM        : vector bien conservative o he toa do local
             *      + TP, TM        : nhiet do
             *      + vectorn       : vector phap tuyen theo he toa do global
             *
            */

            if (DGSchemes::fluxControl::LxF)
                numericalFluxes::LxFFlux(edgeId,flux,fxP,fxM,fyP,fyM,UP,UM,vectorn);
            else if (DGSchemes::fluxControl::central)
                numericalFluxes::centralFlux(flux,fxP,fxM,fyP,fyM,vectorn);
            else if (DGSchemes::fluxControl::Roe)
                numericalFluxes::RoeAverageFlux(flux,fxP,fxM,UP,UM,TP,TM,vectorn);
            else if (DGSchemes::fluxControl::HLL)
                numericalFluxes::HLLFlux(flux,fxP,fxM,UP,UM,TP,TM,vectorn);
        }

        void LxFFlux(int edgeId, std::vector<double>&flux, std::vector<double> &fxP, std::vector<double> &fxM, std::vector<double> &fyP, std::vector<double> &fyM, std::vector<double> &UGP, std::vector<double> &UGM, std::vector<double> &vectorn)
        {
            /* Note: Cac bien su dung trong ham nay la o he toa do Global vi khi su dung LxF flux,
             * khong can xoay ve he toa to Local.
             * Trong ham nay, M la phia ben ngoai, P la phia ben trong phan tu
            */
            double C(LxFConst[edgeId]), nx(vectorn[0]), ny(vectorn[1]);

            for (int i=0; i<4; i++)
            {
                flux[i]=0.5*(nx*(fxP[i]+fxM[i]) + ny*(fyP[i]+fyM[i]) - //chu y dau -
                             C*(UGM[i]-UGP[i]));
            }
        }

        void centralFlux(std::vector<double>&flux, std::vector<double> &fxP, std::vector<double> &fxM, std::vector<double> &fyP, std::vector<double> &fyM, std::vector<double> &vectorn)
        {
            /* Note: Cac bien su dung trong ham nay la o he toa do Global vi khi su dung central flux,
             * khong can xoay ve he toa to Local.
             * Trong ham nay, M la phia ben ngoai, P la phia ben trong phan tu
            */

            double nx(vectorn[0]), ny(vectorn[1]);

            for (int i=0; i<4; i++)
            {
                flux[i]=0.5*(nx*(fxP[i]+fxM[i]) +
                             ny*(fyP[i]+fyM[i]));
            }
        }

        void RoeAverageFlux(std::vector<double>&flux, std::vector<double> &fxPlus, std::vector<double> &fxMinus, std::vector<double> &ULPlus, std::vector<double> &ULMinus, double TPlus, double TMinus, std::vector<double> &vectorn)
        {
            /* Note: Ham nay lay tu sach Nodal Discontinuous Galerkin Methods
             * Trong do quy dinh phia - (Minus, M) la phia ben trong phan tu,
             * phia + (Plus, P) la phia ben ngoai phan tu.
             * Tuy nhien trong code DG quy dinh nguoc lai (+ la ben trong, - la ben ngoai).
             *
             * Vi vay trong ham nay, o phan argument f_Minus la flux cua phia ben
             * ngoai cua phan tu, f_Plus la flux ben trong phan tu.
             * Sau khi truyen vao ham, phai dao nguoc lai de dung voi cong thuc
             * trong sach.
            */
            //Chi tinh flux theo phuong x local vi phuong nay vuong goc voi surface, phuong
            //con lai la tiep tuyen

            //Phai dao nguoc de dung voi cong thuc trong tai lieu
            std::vector<double> ULM(ULPlus), ULP(ULMinus), fxQM(fxPlus), fxQP(fxMinus);

            double uM(ULM[1]/ULM[0]),vM(ULM[2]/ULM[0]),
                    uP(ULP[1]/ULP[0]),vP(ULP[2]/ULP[0]),
                    TM(TPlus),TP(TMinus),
                    rhoM(ULM[0]),rhoP(ULP[0]);

            /*Doan code nay hoan toan tham khao trong sach------------------------------------------------*/
            //Tinh enthalpy: H = totalE + p/rho
            double pM(math::CalcP(TM,rhoM)), pP(math::CalcP(TP,rhoP));

            double HM(material::Cv*TM+0.5*(uM*uM+vM*vM)+pM/rhoM),
                    HP(material::Cv*TP+0.5*(uP*uP+vP*vP)+pP/rhoP);

            //Tinh Roe average variables
            double rhoMs(pow(rhoM,0.5)), rhoPs(pow(rhoP,0.5));

            double
                    rho = rhoMs*rhoPs,
                    u   = (rhoMs*uM + rhoPs*uP)/(rhoMs + rhoPs),
                    v   = (rhoMs*vM + rhoPs*vP)/(rhoMs + rhoPs),
                    H   = (rhoMs*HM + rhoPs*HP)/(rhoMs + rhoPs),

                    c2  = (material::gamma-1)*(H - 0.5*(u*u + v*v)),
                    c   = pow(c2,0.5);

            //Riemann fluxes
            double
                    dW1 = -0.5*rho*(uP-uM)/c + 0.5*(pP-pM)/c2,
                    dW2 = (rhoP-rhoM) - (pP-pM)/c2,
                    dW3 = rho*(vP-vM),
                    dW4 = 0.5*rho*(uP-uM)/c + 0.5*(pP-pM)/c2;

            dW1 = fabs(u-c)*dW1;
            dW2 = fabs(u)*dW2;
            dW3 = fabs(u)*dW3;
            dW4 = fabs(u+c)*dW4;

            //Form Roe fluxes
            std::vector<double> fx(4,0.0);
            for (int i=0; i<4;i++)
            {
                fx[i]=0.5*(fxQP[i]+fxQM[i]);
            }

            fx[0] =fx[0]-(dW1*1       +dW2*1            +dW3*0+dW4*1       )/2;
            fx[1] =fx[1]-(dW1*(u-c)   +dW2*u            +dW3*0+dW4*(u+c)   )/2;
            fx[2] =fx[2]-(dW1*v       +dW2*v            +dW3*1+dW4*v       )/2;
            fx[3] =fx[3]-(dW1*(H-u*c)+dW2*(u*u+v*v)/2+dW3*v+dW4*(H+u*c))/2;

            //rotate back to Cartesian
            flux = fx;
            double nx(vectorn[0]), ny(vectorn[1]);
            flux[1] = nx*fx[1] - ny*fx[2];
            flux[2] = ny*fx[1] + nx*fx[2];
            /*Doan code nay hoan toan tham khao trong sach------------------------------------------------*/
        }

        void HLLFlux(std::vector<double>&flux, std::vector<double> &fxPlus, std::vector<double> &fxMinus, std::vector<double> &ULPlus, std::vector<double> &ULMinus, double TPlus, double TMinus, std::vector<double> &vectorn)
        {
            /* Note: Ham nay lay tu sach Nodal Discontinuous Galerkin Methods
             * Trong do quy dinh phia - (Minus, M) la phia ben trong phan tu,
             * phia + (Plus, P) la phia ben ngoai phan tu.
             * Tuy nhien trong code DG quy dinh nguoc lai (+ la ben trong, - la ben ngoai).
             *
             * Vi vay trong ham nay, o phan argument f_Minus la flux cua phia ben
             * ngoai cua phan tu, f_Plus la flux ben trong phan tu.
             * Sau khi truyen vao ham, phai dao nguoc lai de dung voi cong thuc
             * trong sach.
            */
            //Chi tinh flux theo phuong x local vi phuong nay vuong goc voi surface, phuong
            //con lai la tiep tuyen

            //Phai dao nguoc de dung voi cong thuc trong tai lieu
            std::vector<double> ULM(ULPlus), ULP(ULMinus), fxQM(fxPlus), fxQP(fxMinus);

            double uM(ULM[1]/ULM[0]),vM(ULM[2]/ULM[0]),
                    uP(ULP[1]/ULP[0]),vP(ULP[2]/ULP[0]),
                    TM(TPlus),TP(TMinus),
                    rhoM(ULM[0]),rhoP(ULP[0]);

            /*Doan code nay hoan toan tham khao trong sach------------------------------------------------*/
            //Tinh enthalpy: H = totalE + p/rho
            double pM(math::CalcP(TM,rhoM)), pP(math::CalcP(TP,rhoP));

            double HM(material::Cv*TM+0.5*(uM*uM+vM*vM)+pM/rhoM),
                    HP(material::Cv*TP+0.5*(uP*uP+vP*vP)+pP/rhoP);
            double cM(math::CalcSpeedOfSound(TM)),
                    cP(math::CalcSpeedOfSound(TP));

            //Tinh Roe average variables
            double rhoMs(pow(rhoM,0.5)), rhoPs(pow(rhoP,0.5));

            double
                    //rho = rhoMs*rhoPs,
                    u   = (rhoMs*uM + rhoPs*uP)/(rhoMs + rhoPs),
                    v   = (rhoMs*vM + rhoPs*vP)/(rhoMs + rhoPs),
                    H   = (rhoMs*HM + rhoPs*HP)/(rhoMs + rhoPs),

                    c2  = (material::gamma-1)*(H - 0.5*(u*u + v*v)),
                    c   = pow(c2,0.5);

            //Compute estimate of waves speeds
            double SL = std::min(uM-cM, u-c), SR = std::max(uP+cP, u+c);

            //Riemann fluxes
            /*
            double t1, t2, t3;
            t1 = (std::min(SR,0.0)-std::min(0.0,SL))/(SR-SL);
            t2 = 1-t1;
            t3 = (SR*fabs(SL)-SL*fabs(SR))/(2*(SR-SL));

            //Compute HLL flux
            std::vector<double> fx(4,0.0); //fx means flux
            for (int i=0; i<4;i++)
            {
                fx[i]=t1*fxQP[i] + t2*fxQM[i] - t3*(ULP[i]-ULM[i]);
            }
            */

            std::vector<double> fx(4,0.0); //fx means flux
            if (SL>=0)
            {
                fx=fxQM;
            }
            else if (SR<=0)
            {
                fx=fxQP;
            }
            else
            {
                for (int i=0; i<4;i++)
                {
                    fx[i]=(SR*fxQM[i] - SL*fxQP[i] + SL*SR*(ULP[i]-ULM[i]))/(SR-SL);
                }
            }

            //rotate back to Cartesian
            flux = fx;
            double nx(vectorn[0]), ny(vectorn[1]);
            flux[1] = nx*fx[1] - ny*fx[2];
            flux[2] = ny*fx[1] + nx*fx[2];
            /*Doan code nay hoan toan tham khao trong sach------------------------------------------------*/
        }

        void HLLCFlux(std::vector<double>&flux, std::vector<double> &fxPlus, std::vector<double> &fxMinus, std::vector<double> &ULPlus, std::vector<double> &ULMinus, double TPlus, double TMinus, std::vector<double> &vectorn)
        {
            /* Note: Ham nay chinh sua tu ham HLLFlux lay tu sach Nodal Discontinuous Galerkin Methods
             * Trong do quy dinh phia - (Minus, M) la phia ben trong phan tu,
             * phia + (Plus, P) la phia ben ngoai phan tu.
             * Tuy nhien trong code DG quy dinh nguoc lai (+ la ben trong, - la ben ngoai).
             *
             * Vi vay trong ham nay, o phan argument f_Minus la flux cua phia ben
             * ngoai cua phan tu, f_Plus la flux ben trong phan tu.
             * Sau khi truyen vao ham, phai dao nguoc lai de dung voi cong thuc
             * trong sach.
            */
            //Chi tinh flux theo phuong x local vi phuong nay vuong goc voi surface, phuong
            //con lai la tiep tuyen

            //Phai dao nguoc de dung voi cong thuc trong tai lieu
            std::vector<double> ULM(ULPlus), ULP(ULMinus), fxQM(fxPlus), fxQP(fxMinus);

            double uM(ULM[1]/ULM[0]),vM(ULM[2]/ULM[0]),
                    uP(ULP[1]/ULP[0]),vP(ULP[2]/ULP[0]),
                    TM(TPlus),TP(TMinus),
                    rhoM(ULM[0]),rhoP(ULP[0]);

            /*Doan code nay hoan toan tham khao trong sach------------------------------------------------*/
            //Tinh enthalpy: H = totalE + p/rho
            double pM(math::CalcP(TM,rhoM)), pP(math::CalcP(TP,rhoM));

            double HM(material::Cv*TM+0.5*(uM*uM+vM*vM)+pM/rhoM),
                    HP(material::Cv*TP+0.5*(uP*uP+vP*vP)+pP/rhoM);
            double cM(math::CalcSpeedOfSound(TM)),
                    cP(math::CalcSpeedOfSound(TP));

            //Tinh Roe average variables
            double rhoMs(pow(rhoM,0.5)), rhoPs(pow(rhoP,0.5));

            double
                    //rho = rhoMs*rhoPs,
                    u   = (rhoMs*uM + rhoPs*uP)/(rhoMs + rhoPs),
                    v   = (rhoMs*vM + rhoPs*vP)/(rhoMs + rhoPs),
                    H   = (rhoMs*HM + rhoPs*HP)/(rhoMs + rhoPs),

                    c2  = (material::gamma-1)*(H - 0.5*(u*u + v*v)),
                    c   = pow(c2,0.5);

            //Compute estimate of waves speeds
            double SL = std::min(uM-cM, u-c), SR = std::max(uP+cP, u+c);

            //Riemann fluxes
            /*
            double t1, t2, t3;
            t1 = (std::min(SR,0.0)-std::min(0.0,SL))/(SR-SL);
            t2 = 1-t1;
            t3 = (SR*fabs(SL)-SL*fabs(SR))/(2*(SR-SL));

            //Compute HLL flux
            std::vector<double> fx(4,0.0); //fx means flux
            for (int i=0; i<4;i++)
            {
                fx[i]=t1*fxQP[i] + t2*fxQM[i] - t3*(ULP[i]-ULM[i]);
            }
            */

            std::vector<double> fx(4,0.0); //fx means flux
            if (SL>=0)
            {
                fx=fxQM;
            }
            else if (SR<=0)
            {
                fx=fxQP;
            }
            else
            {
                for (int i=0; i<4;i++)
                {
                    fx[i]=(SR*fxQM[i] - SL*fxQP[i] + SL*SR*(ULP[i]-ULM[i]))/(SR-SL);
                }
            }

            //rotate back to Cartesian
            flux = fx;
            double nx(vectorn[0]), ny(vectorn[1]);
            flux[1] = nx*fx[1] - ny*fx[2];
            flux[2] = ny*fx[1] + nx*fx[2];
            /*Doan code nay hoan toan tham khao trong sach------------------------------------------------*/
        }

        //Dung de tinh flux cho phuong trinh phu, tuong duong central flux
        double auxFlux(double MinusVal, double PlusVar, double vectorComp)
        {
            /*use central numerical flux*/
            double flux(0.5*(MinusVal + PlusVar)*vectorComp);
            return flux;
        }

        double constantC(double uMagP, double uMagM, double aP, double aM)
        {
            std::vector<double> CArray(2, 0.0);
            double C(0.0);
            CArray[0] = uMagP + aP;
            CArray[1] = uMagM + aM;
            C = *std::max_element(CArray.begin(), CArray.end());  //find max value of vector
            return C;
        }

        void findMaxLxFConstantOnEdge(int iedge, int masterCell)
        {
            /* Ham tinh he so C tren 1 edge trong truong hop Flux type la LxF,
             * cell id input vao ham la master cell cua edge
            */
            double TMaster, TSlave;
            std::vector<double> UMaster(4), USlave(4);
            std::vector<double> CArray(mathVar::nGauss1D+1);
            for (int nG = 0; nG <= mathVar::nGauss1D; nG++)
            {
                std::tie(TMaster,TSlave)=auxUlti::getTAtInterfaces(iedge,masterCell,nG);
                for (int i = 0; i < 4; i++)
                {
                    //std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
                    std::tie(UMaster[i], USlave[i]) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, i + 1);
                }
                double uMagMaster(sqrt(pow(UMaster[1]/UMaster[0],2)+pow(UMaster[2]/UMaster[0],2))),
                        uMagSlave(sqrt(pow(USlave[1]/USlave[0],2)+pow(USlave[2]/USlave[0],2))),
                        aMaster(math::CalcSpeedOfSound(TMaster)),
                        aSlave(math::CalcSpeedOfSound(TSlave));

                CArray[nG]=(math::numericalFluxes::constantC(uMagMaster,uMagSlave,aMaster,aSlave));
            }
            LxFConst[iedge]=*max_element(CArray.begin(), CArray.end());
            //-----------------------------------------------
        }
    }

    namespace inviscidTerms
    {
        std::tuple<double, double, double, double> calcInvisTermsFromPriVars(double rhoVal, double uVal, double umVal, double vVal, double vmVal, double totalE, double pVal, int dir)
        {
            double term1(0.0), term2(0.0), term3(0.0), term4(0.0);

            if (dir==1)  //Ox direction
            {
                term1 = rhoVal * uVal;
                term2 = rhoVal * uVal*umVal + pVal;
                term3 = rhoVal * uVal*vmVal;
                term4 = rhoVal*totalE*umVal + pVal*uVal;
            }
            else if (dir==2)  //Oy direction
            {
                term1 = rhoVal * vVal;
                term2 = rhoVal * umVal*vVal;
                term3 = rhoVal * vVal*vmVal + pVal;
                term4 = rhoVal*totalE*vmVal + pVal*vVal;
            }
            return std::make_tuple(term1, term2, term3, term4);
        }
    }

    namespace viscousTerms
    {
        std::vector<std::vector<double>> calcStressTensorAndHeatFlux(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
        {
            /*Output matrix has form:
            [tauXx		tauXy		Qx]
            [tauYx		tauYy		Qy]
            */
            std::vector<std::vector<double>> OutputMatrix(2, std::vector<double>(3, 0.0));

            double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), dEx(0.0), dEy(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0), rhoEVal(0.0);
            double drhox(0.0), drhoy(0.0),
                drhoux(0.0), drhouy(0.0),
                drhovx(0.0), drhovy(0.0),
                drhoEx(0.0), drhoEy(0.0),
                dTx(0.0), dTy(0.0);
            double uVal(0.0), vVal(0.0), Qx(0.0), Qy(0.0), k(0.0), uMag(0.0);
            int index(0);

            rhoVal = U[0];
            rhouVal = U[1];
            rhovVal = U[2];

            rhoEVal = U[3];

            uVal = rhouVal / rhoVal;
            vVal = rhovVal / rhoVal;
            uMag=pow((uVal*uVal+vVal*vVal),0.5);

            drhox = dUx[0];
            drhoy = dUy[0];

            drhoux = dUx[1];
            drhouy = dUy[1];

            drhovx = dUx[2];
            drhovy = dUy[2];

            drhoEx = dUx[3];
            drhoEy = dUy[3];

            dux = math::calcRhouvEDeriv(drhoux, drhox, rhouVal, rhoVal);
            duy = math::calcRhouvEDeriv(drhouy, drhoy, rhouVal, rhoVal);

            dvx = math::calcRhouvEDeriv(drhovx, drhox, rhovVal, rhoVal);
            dvy = math::calcRhouvEDeriv(drhovy, drhoy, rhovVal, rhoVal);

            dEx = math::calcRhouvEDeriv(drhoEx, drhox, rhoEVal, rhoVal);
            dEy = math::calcRhouvEDeriv(drhoEy, drhoy, rhoEVal, rhoVal);

            /*calculate stresses*/
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    index = 10 * (i + 1) + (j + 1);
                    if (index == 11)  //x normal stress (tau_xx)
                    {
                        OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, dux, dvy);
                    }
                    else if (index == 22)  //y normal stress (tau_yy)
                    {
                        OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, dvy, dux);
                    }
                    else  //shear stress (tau_xy)
                    {
                        OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, duy, dvx);
                    }
                }
            }

            /*calculate heat flux*/
            dTx = math::calcTDeriv(dEx, dux, dvx, uVal, vVal);
            dTy = math::calcTDeriv(dEy, duy, dvy, uVal, vVal);
            //k = math::calcThermalConductivity(muVal);
            //dTx=math::calcDivTFromAuxVariable(drhoEx,drhox,rhoEVal,rhoVal);
            //dTy=math::calcDivTFromAuxVariable(drhoEy,drhoy,rhoEVal,rhoVal);
            k = material::Cp / material::Pr;
            std::tie(Qx, Qy) = math::viscousTerms::calcHeatFluxTerms(dTx, dTy, k);

            OutputMatrix[0][2] = Qx;
            OutputMatrix[1][2] = Qy;
            return OutputMatrix;
        }

        double calcStressComponent(int index, double fstDeriv, double sndDeriv)
        {
            /*Formula of stress
            tau_xx = -mu((4/3)*du/dx - (2/3)*dv/dy)
            tau_yy = -mu((4/3)*dv/dy - (2/3)*du/dx)
            tau_xy = -mu(du/dy + dv/dx)

            input variables:
            index			fstDeriv		sndDeriv
            11 (tau_xx)		du/dx			dv/dy
            22 (tau_yy)		dv/dy			du/dx
            12 (tau_xy)		du/dy			dv/dx
            */
            double tau(0.0);
            if (index == 11 || index == 22)
            {
                tau = -((4.0 / 3.0)*fstDeriv - (2.0 / 3.0)*sndDeriv);
            }
            else if (index == 12 || index == 21)
            {
                tau = -(fstDeriv + sndDeriv);
            }
            return tau;
        }

        std::tuple<double, double> calcHeatFluxTerms(double dTx, double dTy, double k)
        {
            double Qx(0.0), Qy(0.0);
            Qx = -k * dTx;
            Qy = -k * dTy;
            return std::make_tuple(Qx, Qy);
        }

        std::tuple<double, double, double, double> calcViscousTermsFromStressHeatFluxMatrix(std::vector< std::vector<double> > &StressHeatFlux, double uVal, double vVal, int dir)
        {
            /*StressHeatFlux is 2D array contents stress and heat flux component, has the following form:
            [tauXx		tauXy		Qx]
            [tauYx		tauYy		Qy]*/

            double tauXy(0.0), tauXx(0.0), tauYy(0.0), Qx(0.0), Qy(0.0);
            double viscTerm1(0.0), viscTerm2(0.0), viscTerm3(0.0), viscTerm4(0.0);

            tauXx = StressHeatFlux[0][0];
            tauYy = StressHeatFlux[1][1];
            tauXy = StressHeatFlux[0][1];
            Qx = StressHeatFlux[0][2];
            Qy = StressHeatFlux[1][2];

            if (dir==1)
            {
                /*1. Ox direction*/
                viscTerm1=0.0;
                viscTerm2 = tauXx;
                viscTerm3 = tauXy;
                viscTerm4 = tauXx * uVal + tauXy * vVal + Qx;
            }
            else if (dir==2)
            {
                /*2. Oy direction*/
                viscTerm1=0.0;
                viscTerm2 = tauXy;
                viscTerm3 = tauYy;
                viscTerm4 = tauXy * uVal + tauYy * vVal + Qy;
            }
            return std::make_tuple(viscTerm1, viscTerm2, viscTerm3, viscTerm4);
        }
    }//end of namespace viscousTerms

    /*Function computes value of conservative variables at abitrary point with applying limiter
            valType:
            1: rho
            2: rhou
            3: rhov
            4: rhoE*/
    double calcConsvVarWthLimiter(int element, double a, double b, int valType)
    {
        double out(0.0);
        std::vector<double> Value(mathVar::orderElem + 1, 0.0);
        Value = auxUlti::getElementConserValuesOfOrder(element, valType);

        //Compute value at point (a, b) without limiter
        math::basisFc(a, b, auxUlti::checkType(element));

        if (valType == 1)
        {
            for (int order = 1; order <= mathVar::orderElem; order++)
            {
                out += Value[order] * mathVar::B[order] * theta2Arr[element] * theta1Arr[element];
            }
        }
        else
        {
            for (int order = 1; order <= mathVar::orderElem; order++)
            {
                out += Value[order] * mathVar::B[order] * theta2Arr[element];
            }
        }
        out += Value[0];

        return out;
    }

    double calcResidualFromResidualOfOrder(int element, double a, double b, int valType)
    {
        double out(0.0);
        std::vector<double> Value(mathVar::orderElem + 1, 0.0);
        Value = auxUlti::getResidualValuesOfOrder(element, valType);

        //Compute value at point (a, b) without limiter
        math::basisFc(a, b, auxUlti::checkType(element));

        for (int order = 1; order <= mathVar::orderElem; order++)
        {
            out += Value[order] * mathVar::B[order];
        }
        out += Value[0];

        return fabs(out);
    }

    namespace geometricOp
    {
        double calcDistanceFromCenterToEdge(int element, int edge)
        {
            double //Coordinates of centroid
                    xC=meshVar::geoCenter[element][0],
                    yC=meshVar::geoCenter[element][1],

                    //Coordinates of 2 vertices of edge
                    x1=meshVar::Points[meshVar::inpoed[edge][0]][0],
                    y1=meshVar::Points[meshVar::inpoed[edge][0]][1],
                    x2=meshVar::Points[meshVar::inpoed[edge][1]][0],
                    y2=meshVar::Points[meshVar::inpoed[edge][1]][1];
            //Function of line: y=ax+b <=> ax-y+b=0
            double a((y1-y2)/(x1-x2)),
                    b=(-((y1-y2)/(x1-x2))*x2+y2);

            //Distance from centroid to edge
            return (fabs(a*xC-yC+b)/pow((a*a+1),0.5));
        }

        std::tuple<double, double> calcGeoCenter(std::vector<double> &xCoor, std::vector<double> &yCoor, int type)
        {
            double x(0.0), y(0.0), xCG(0.0), yCG(0.0);
            for (int i = 0; i < type; i++)
            {
                x = xCoor[i];
                y = yCoor[i];
                xCG += x;
                yCG += y;
            }

            xCG = xCG / type;
            yCG = yCG / type;
            return std::make_tuple(xCG, yCG);
        }

        double calcPolygonArea(std::vector<double> &xCoor, std::vector<double> &yCoor, int type)
        {
            double x1(xCoor[0]), x2(xCoor[1]), x3(xCoor[2]), x4(0.0), y1(yCoor[0]), y2(yCoor[1]), y3(yCoor[2]), y4(0.0),Area(0.0);
            if (type == 3)
            {
                Area = fabs(0.5*(x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3));
            }
            else if (type == 4)
            {
                x4 = xCoor[3];
                y4 = yCoor[3];

                Area = fabs(0.5*(x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y4 - x4 * y3 + x4 * y1 - x1 * y4));
            }
            return Area;
        }

        std::tuple<double, double> calcQuadCentroid(int element, double xCG, double yCG, double area)
        {
            std::vector<double> xSubTriCoor(3, 0.0), ySubTriCoor(3, 0.0);
            std::vector<double> xCGSubTri(4, 0.0), yCGSubTri(4, 0.0), subTriArea(4,0.0);
            double xC(0.0), yC(0.0);
            //1. point 0, 1
            std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 0);
            std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 1);
            xSubTriCoor[2] = xCG;
            ySubTriCoor[2] = yCG;
            std::tie(xCGSubTri[0], yCGSubTri[0]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
            subTriArea[0] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

            //2. point 1, 2
            std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 1);
            std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 2);
            xSubTriCoor[2] = xCG;
            ySubTriCoor[2] = yCG;
            std::tie(xCGSubTri[1], yCGSubTri[1]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
            subTriArea[1] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

            //3. point 2, 3
            std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 2);
            std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 3);
            xSubTriCoor[2] = xCG;
            ySubTriCoor[2] = yCG;
            std::tie(xCGSubTri[2], yCGSubTri[2]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
            subTriArea[2] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

            //4. point 3, 0
            std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 3);
            std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 0);
            xSubTriCoor[2] = xCG;
            ySubTriCoor[2] = yCG;
            std::tie(xCGSubTri[3], yCGSubTri[3]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
            subTriArea[3] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

            for (int i = 0; i < 4; i++)
            {
                xC += subTriArea[i] * xCGSubTri[i];
                yC += subTriArea[i] * yCGSubTri[i];
            }
            xC = xC / area;
            yC = yC / area;
            return std::make_tuple(xC, yC);
        }

        double calLocalCellSize(int element, double elementArea)
        {
            int elemType(auxUlti::checkType(element)), edgeId(0), neighborElemId(0);
            double deltaXe(0.0), deltaYe(0.0), Lx(0.0), Ly(0.0), Lxy(0.0), deltaXc(0.0), deltaYc(0.0);

            for (int e = 0; e < elemType; e++)
            {
                edgeId = meshVar::inedel[element][e];
                neighborElemId = auxUlti::getNeighborElement(element, edgeId);
                std::tie(deltaXe, deltaYe) = math::geometricOp::calEdgeMetric(edgeId, element);
                if (neighborElemId >= 0)
                {
                    std::tie(deltaXc, deltaYc) = math::geometricOp::calDifferenceOfElementsCellCenters(element, neighborElemId);
                }
                else
                {
                    deltaXc = 0.0;
                    deltaYc = 0.0;
                }

                Lx += deltaXc * fabs(deltaXc)*deltaYe;
                Ly += deltaYc * fabs(deltaYc)*deltaXe;
            }
            Lx *= (1.0 / elementArea);
            Ly *= (-1.0 / elementArea);
            Lxy = sqrt(Lx*Lx + Ly*Ly);
            return Lxy;
        }

        std::tuple<double, double> calEdgeMetric(int edgeId, int elementId)
        {
            int point1(meshVar::inpoed[edgeId][0]), point2(meshVar::inpoed[edgeId][1]);
            double deltaX(0.0), deltaY(0.0);
            if (auxUlti::findVertexOrder(point1, elementId) > auxUlti::findVertexOrder(point2, elementId))
            {
                deltaX = meshVar::Points[point1][0] - meshVar::Points[point2][0];
                deltaY = meshVar::Points[point1][1] - meshVar::Points[point2][1];
            }
            else
            {
                deltaX = meshVar::Points[point2][0] - meshVar::Points[point1][0];
                deltaY = meshVar::Points[point2][1] - meshVar::Points[point1][1];
            }
            return std::make_tuple(deltaX, deltaY);
        }

        std::tuple<double, double> calDifferenceOfElementsCellCenters(int elem1, int elem2)
        {
            double deltaXc(0.0), deltaYc(0.0);
            //elem2 is neighbor
            deltaXc = (meshVar::geoCenter[elem2][0] - meshVar::geoCenter[elem1][0]);
            deltaYc = (meshVar::geoCenter[elem2][1] - meshVar::geoCenter[elem1][1]);
            return std::make_tuple(deltaXc, deltaYc);
        }

        double calDistBetween2Points(double xPt1, double yPt1, double xPt2, double yPt2) {
            return (sqrt(pow(xPt1 - xPt2, 2) + pow(yPt1 - yPt2, 2)));
        }

        std::tuple<double, double> calcUnitVectorOf2Points(double xPt1, double yPt1, double xPt2, double yPt2)
        {
            double length(sqrt(pow(xPt2 - xPt1, 2) + pow(yPt2 - yPt1, 2)));
            return std::make_tuple
                    (
                        (xPt2 - xPt1)/length,
                        (yPt2 - yPt1)/length
                    );
        }

        double calPolygonPerimeter(std::vector<double> &xCoor, std::vector<double> &yCoor, int numOfEdge) {
            double perimeter(0.0);
            for (int iedge = 0; iedge < numOfEdge - 1; iedge++)
            {
                perimeter += math::geometricOp::calDistBetween2Points(xCoor[iedge], yCoor[iedge], xCoor[iedge + 1], yCoor[iedge + 1]);
            }
            perimeter += math::geometricOp::calDistBetween2Points(xCoor[numOfEdge -1], yCoor[numOfEdge - 1], xCoor[0], yCoor[0]);
            return perimeter;
        }

        std::tuple<double, double> calROfInscribedCircleOfTriElement(std::vector<double> &xCoor, std::vector<double> &yCoor) {
            double rIn(0.0), cellArea(0.0);
            cellArea = math::geometricOp::calcPolygonArea(xCoor, yCoor, 3);
            rIn = 2 * cellArea / math::geometricOp::calPolygonPerimeter(xCoor, yCoor, 3);
            return std::make_tuple(rIn, cellArea);
        }

        std::tuple<double, double> calSizeOfQuadElement(std::vector<double> &xCoor, std::vector<double> &yCoor) {
            int elemType(4);
            double size(0.0), size1(0.0), size2(0.0), cellArea(0.0);
            cellArea = math::geometricOp::calcPolygonArea(xCoor, yCoor, elemType);
            size1 = math::geometricOp::calDistBetween2Points(xCoor[0], yCoor[0], xCoor[2], yCoor[2]);
            size2 = math::geometricOp::calDistBetween2Points(xCoor[1], yCoor[1], xCoor[3], yCoor[3]);
            if (size1>size2)
            {
                size=size2;
            }
            else
            {
                size=size1;
            }
            return std::make_tuple(size, cellArea);
        }

        std::tuple<double, double> findNormalProjectionOfPointToEdge(double xA, double yA, double xB, double yB, double xC, double yC)
        {
            //Project point A to edge BC
            double xN, yN;

            if (xB==xC && yB!=yC)
            {
                xN = xB;
                yN = yA;
            }
            else if (xB!=xC && yB==yC)
            {
                xN = xA;
                yN = yB;
            }
            else if (xB!=xC && yB!=yC)
            {
                double a((yB-yC)/(xB-xC)),
                        b=(-((yB-yC)/(xB-xC))*xB+yB);
                //Phuong trinh duong thang la y=ax+b
                //Phuong trinh duong thang vuong goc la y=-x/a+bN
                double bN(yA+xA/a);
                //Toa do hinh chieu vuong goc cua A len BC
                xN = (bN-b)/(a+1/a);
                yN = a*xN+b;
            }
            else
            {
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, "Mesh error: 2 vertices of 1 edges are coincident.");
            }

            return std::make_tuple(xN, yN);
        }

        std::vector<int> sortVerticesCCW(std::vector<int> verticesId, std::vector<double> angle)
        {
            /* Ham sap xep cac vertex cua 1 cell theo chieu counter clockwise.
             * Argument la verticesId cu, goc tao boi vertex, tam cua hinh poly va
             * truc Ox cua tung vertex. Ham se sap xep cac vertex theo chieu goc
             * tang dan
            */
            int nV(static_cast<int>(verticesId.size())), id;
            std::vector<int> newOrder(nV,0),
                    sortVector(nV,0);
            double m;

            for (int i=0; i<nV; i++)
            {
                m=*max_element(angle.begin(), angle.end());
                id=math::findIndex_double(m,angle);
                sortVector[i]=id;
                angle[id]=370.0;
            }

            for (int i=0; i<nV; i++)
            {
                id=sortVector[i];
                newOrder[i]=verticesId[id];
            }
            return newOrder;
        }

        double calcAngleOfPoint(double xOrig,double yOrig,double xPoint,double yPoint)
        {
            /* Ham tinh goc tao boi 1 dinh cua poly, tam poly va chieu + truc Ox
            */
            double A[2]={fabs(1.5*(xOrig+0.1)), yOrig},
                    O[2]={xOrig,yOrig},
                    B[2]={xPoint,yPoint};

            double OA[2]={A[0]-O[0],A[1]-O[1]},
                    OB[2]={B[0]-O[0],B[1]-O[1]};

            double lOA(sqrt(OA[0]*OA[0]+OA[1]*OA[1])),
                    lOB(sqrt(OB[0]*OB[0]+OB[1]*OB[1])),

                    dotProduct=OA[0]*OB[0]+OA[1]*OB[1],
                    det=OA[0]*OB[1]-OA[1]*OB[0];

            double cosx=dotProduct/(lOA*lOB), rad, angle;
            rad=acos(cosx);
            angle=rad*180/3.1416;

            if (det>0) return angle;
            else return (360-angle);
        }

        void findGaussPtsOnWallParameters()
        {
            int globleEdge(0), element(0);
            double a, b, xC, yC, xGauss, yGauss;
            for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
            {
                globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
                //int bcType(auxUlti::checkBCTypeOfEdge(globleEdge));
                std::tie(element,std::ignore)=auxUlti::getMasterServantOfEdge(globleEdge);
                for (int nG=0; nG<=mathVar::nGauss1D; nG++)
                {
                    //Get Gauss point's coordinates
                    std::tie(a, b) = auxUlti::getGaussSurfCoor(globleEdge, element, nG);

                    //Map from standard to real coor sys
                    std::tie(xGauss, yGauss) = math::directMapping(element, a, b);

                    meshVar::GaussPtsOnBCEdge_x[ilocalEdge][nG] = xGauss;
                    meshVar::GaussPtsOnBCEdge_y[ilocalEdge][nG] = yGauss;

                    //Get centroid
                    xC=meshVar::geoCenter[element][0];
                    yC=meshVar::geoCenter[element][1];

                    //Calc distance from centroid to Gauss point
                    meshVar::distanceFromGaussPtsToCentroid[ilocalEdge][nG] = math::geometricOp::calDistBetween2Points(xC,yC,xGauss,yGauss);

                    //Find unit vector (Centroid)---------->(Gauss on edge)
                    std::tie(meshVar::GaussPtsOnBCEdge_unitVector_x[ilocalEdge][nG],
                             meshVar::GaussPtsOnBCEdge_unitVector_y[ilocalEdge][nG])
                            =math::geometricOp::calcUnitVectorOf2Points(xC,yC,xGauss,yGauss);
                }
            }
        }
    }

    namespace residualManipulation
    {
        /*
        void calcNormResidual(double rhoRes, double rhouRes, double rhovRes, double rhoERes)
        {
            systemVar::rhoResNormVector[systemVar::iterCount - 1] = rhoRes;
            systemVar::rhouResNormVector[systemVar::iterCount - 1] = rhouRes;
            systemVar::rhovResNormVector[systemVar::iterCount - 1] = rhovRes;
            systemVar::rhoEResNormVector[systemVar::iterCount - 1] = rhoERes;

            systemVar::rhoResNorm = *std::max_element(systemVar::rhoResNormVector.begin(), systemVar::rhoResNormVector.end());  //find min value of vector
            systemVar::rhouResNorm = *std::max_element(systemVar::rhouResNormVector.begin(), systemVar::rhouResNormVector.end());
            systemVar::rhovResNorm = *std::max_element(systemVar::rhovResNormVector.begin(), systemVar::rhovResNormVector.end());
            systemVar::rhoEResNorm = *std::max_element(systemVar::rhoEResNormVector.begin(), systemVar::rhoEResNormVector.end());
        }
        */
    }

    double calcMaxT(int element)
    {
        std::vector<double> vectorT(2 * (mathVar::nGauss2D + 1) * (mathVar::nGauss2D + 1), 0.0);
        double aG(0.0), bG(0.0), aGL(0.0), bGL(0.0), min(0.0);
        int index(0);
        for (int na = 0; na <= mathVar::nGauss2D; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss2D; nb++)
            {
                int nanb(calcArrId(na,nb,mathVar::nGauss2D+1));
                std::tie(aG,bG)=auxUlti::getGaussCoor(na,nb);

                aGL = mathVar::GaussLobattoPts[nanb][0];
                bGL = mathVar::GaussLobattoPts[nanb][1];

                vectorT[index] = math::pointValue(element, aG, bGL, 6, 1);
                index++;
                vectorT[index] = math::pointValue(element, aGL, bG, 6, 1);
                index++;
            }
        }
        min = *std::max_element(vectorT.begin(), vectorT.end());  //find min value of vector
        return min;
    }

    int findIndex(int number, std::vector<int> InArray)
    {
        int index(0), size(static_cast<int>(InArray.size()));
        for (int i = 0; i < size; i++)
        {
            if (number == InArray[i])
            {
                index = i;
                break;
            }
            else
            {
                index = -1;
            }
        }
        return index;
    }

    int findIndex_double(double number, std::vector<double> InArray)
    {
        int index(0), size(static_cast<int>(InArray.size()));
        for (int i = 0; i < size; i++)
        {
            if (fabs(InArray[i]-number)<1e-10)
            {
                index = i;
                break;
            }
            else
            {
                index = -1;
            }
        }
        return index;
    }

    std::vector<double> pointUVars(int element, double a, double b)
    {
        std::vector<double> U(4, 0.0);

        math::basisFc(a, b, auxUlti::checkType(element));
        for (int iorder = 1; iorder <= mathVar::orderElem; iorder++)
        {
            U[0] += rho[element][iorder]* mathVar::B[iorder] * theta2Arr[element] * theta1Arr[element];
            U[1] += rhou[element][iorder]* mathVar::B[iorder] * theta2Arr[element];
            U[2] += rhov[element][iorder]* mathVar::B[iorder] * theta2Arr[element];
            U[3] += rhoE[element][iorder]* mathVar::B[iorder] * theta2Arr[element];
        }
        U[0] += rho[element][0];
        U[1] += rhou[element][0];
        U[2] += rhov[element][0];
        U[3] += rhoE[element][0];
        return U;
    }

    std::vector<double> pointSVars(int edge, int element, double a, double b, int dir, int option)
    {
        //option dung trong truong hop BRScheme la BR2:
        //- option=1: volume
        //- option=2: surface
        std::vector<double> dU(4, 0.0);

        math::basisFc(a, b, auxUlti::checkType(element));
        //Limiter
        double theta(1.0);
        if (flowProperties::massDiffusion)
        {
            theta=theta1Arr[element]*theta2Arr[element];
        }

        if (systemVar::auxVariables==1)
        {
            if (dir==1)
            {
                for (int iorder = 1; iorder <= mathVar::orderElem; iorder++)
                {
                    dU[0] += BR1Vars::rhoX[element][iorder]* mathVar::B[iorder]*theta;
                    dU[1] += BR1Vars::rhouX[element][iorder]* mathVar::B[iorder]*theta;
                    dU[2] += BR1Vars::rhovX[element][iorder]* mathVar::B[iorder]*theta;
                    dU[3] += BR1Vars::rhoEX[element][iorder]* mathVar::B[iorder]*theta;
                }
                dU[0] += BR1Vars::rhoX[element][0];
                dU[1] += BR1Vars::rhouX[element][0];
                dU[2] += BR1Vars::rhovX[element][0];
                dU[3] += BR1Vars::rhoEX[element][0];
            }
            else {
                for (int iorder = 1; iorder <= mathVar::orderElem; iorder++)
                {
                    dU[0] += BR1Vars::rhoY[element][iorder]* mathVar::B[iorder]*theta;
                    dU[1] += BR1Vars::rhouY[element][iorder]* mathVar::B[iorder]*theta;
                    dU[2] += BR1Vars::rhovY[element][iorder]* mathVar::B[iorder]*theta;
                    dU[3] += BR1Vars::rhoEY[element][iorder]* mathVar::B[iorder]*theta;
                }
                dU[0] += BR1Vars::rhoY[element][0];
                dU[1] += BR1Vars::rhouY[element][0];
                dU[2] += BR1Vars::rhovY[element][0];
                dU[3] += BR1Vars::rhoEY[element][0];
            }
        }
        else if (systemVar::auxVariables==2)
        {
            if (option==1) //volume
            {
                if (dir==1)
                {
                    for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                    {
                        dU[0] += BR2Vars::rhoXVol[element][iorder]* mathVar::B[iorder];
                        dU[1] += BR2Vars::rhouXVol[element][iorder]* mathVar::B[iorder];
                        dU[2] += BR2Vars::rhovXVol[element][iorder]* mathVar::B[iorder];
                        dU[3] += BR2Vars::rhoEXVol[element][iorder]* mathVar::B[iorder];
                    }
                }
                else {
                    for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                    {
                        dU[0] += BR2Vars::rhoYVol[element][iorder]* mathVar::B[iorder];
                        dU[1] += BR2Vars::rhouYVol[element][iorder]* mathVar::B[iorder];
                        dU[2] += BR2Vars::rhovYVol[element][iorder]* mathVar::B[iorder];
                        dU[3] += BR2Vars::rhoEYVol[element][iorder]* mathVar::B[iorder];
                    }
                }
            }
            else if (option==2) { //surface
                if (auxUlti::checkMaster(element,edge)) //element is master of edge
                {
                    if (dir==1)
                    {
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            dU[0] += BR2Vars::rhoXSurMaster[edge][iorder]* mathVar::B[iorder];
                            dU[1] += BR2Vars::rhouXSurMaster[edge][iorder]* mathVar::B[iorder];
                            dU[2] += BR2Vars::rhovXSurMaster[edge][iorder]* mathVar::B[iorder];
                            dU[3] += BR2Vars::rhoEXSurMaster[edge][iorder]* mathVar::B[iorder];
                        }
                    }
                    else {
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            dU[0] += BR2Vars::rhoYSurMaster[edge][iorder]* mathVar::B[iorder];
                            dU[1] += BR2Vars::rhouYSurMaster[edge][iorder]* mathVar::B[iorder];
                            dU[2] += BR2Vars::rhovYSurMaster[edge][iorder]* mathVar::B[iorder];
                            dU[3] += BR2Vars::rhoEYSurMaster[edge][iorder]* mathVar::B[iorder];
                        }
                    }
                }
                else {
                    if (dir==1)
                    {
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            dU[0] += BR2Vars::rhoXSurSlave[edge][iorder]* mathVar::B[iorder];
                            dU[1] += BR2Vars::rhouXSurSlave[edge][iorder]* mathVar::B[iorder];
                            dU[2] += BR2Vars::rhovXSurSlave[edge][iorder]* mathVar::B[iorder];
                            dU[3] += BR2Vars::rhoEXSurSlave[edge][iorder]* mathVar::B[iorder];
                        }
                    }
                    else {
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            dU[0] += BR2Vars::rhoYSurSlave[edge][iorder]* mathVar::B[iorder];
                            dU[1] += BR2Vars::rhouYSurSlave[edge][iorder]* mathVar::B[iorder];
                            dU[2] += BR2Vars::rhovYSurSlave[edge][iorder]* mathVar::B[iorder];
                            dU[3] += BR2Vars::rhoEYSurSlave[edge][iorder]* mathVar::B[iorder];
                        }
                    }
                }
            }
        }

        return dU;
    }

    std::tuple<double, double> rotateToLocalCoordinateSystem(double xComp, double yComp, double nx, double ny)
    {
        /* localXComp tuong duong Ox
         * localYComp tuong duong Oy
        */
        double localXComp(nx*xComp+ny*yComp);
        double localYComp(-ny*xComp+nx*yComp);
        return std::make_tuple(localXComp,localYComp);
    }

    std::tuple<double, double> rotateToGlobalCoordinateSystem(double localXComp, double localYComp, double nx, double ny)
    {
        /* normalComp tuong duong Ox
         * tangentialComp tuong duong Oy
        */
        double xComp(nx*localXComp-ny*localYComp);
        double yComp(ny*localXComp+nx*localYComp);
        return std::make_tuple(xComp,yComp);
    }

    namespace solvePolynomialsEq {
    double NewtonRaphson(std::vector<double> &power, std::vector<double> &coefs, double initialValue)
        {
            //tang initial value len
            double initialValue_org(initialValue);
            initialValue=initialValue*5.0;

            int polySize(static_cast<int>(power.size())), maxIter(50);
            std::vector<double> power_deriv(polySize, 0.0),
                    coefs_deriv(polySize, 0.0);
            double error(1), convergence_cri(1e-8), output(0.0), fValue(0.0), derivfValue(0.0);

            //find 1st derivative of input polynomial
            for (int polyOrder = 0; polyOrder < polySize; polyOrder++) {
                power_deriv[polyOrder]=power[polyOrder] - 1.0;
                coefs_deriv[polyOrder]=power[polyOrder]*coefs[polyOrder];
            }

            //solve equation
            int iter(0);
            while (error > convergence_cri && iter<maxIter) {
                fValue = 0.0;
                derivfValue = 0.0;
                for (int polyOrder = 0; polyOrder < polySize; polyOrder++) {
                    fValue+=pow(initialValue,power[polyOrder])*coefs[polyOrder];
                    derivfValue+=pow(initialValue,power_deriv[polyOrder])*coefs_deriv[polyOrder];
                }
                output = fabs(initialValue - fValue/derivfValue);
                error = fabs(output - initialValue)/initialValue;
                initialValue = output;
                iter++;
            }

            if (error<convergence_cri)
            {
                return output;
            }
            else
            {
                return initialValue_org;
            }
        }

        double Bisection(std::vector<double> &power, std::vector<double> &coefs, double initialValue)
        {
            const double range(0.3);
            double upper(initialValue*(1+range)), lower((initialValue*(1-range))), fUpper(0.0), fLower(0.0);
            double error(1), convergence_cri(1e-6), output(0.0), fOutput(0.0), outputOld(0.0);
            //solve equation
            fUpper=math::solvePolynomialsEq::subValToPolynomial(power,coefs,upper);
            fLower=math::solvePolynomialsEq::subValToPolynomial(power,coefs,lower);
            if (fUpper*fLower>0)
            {
                std::cout<<"Failed to solve polynomial equation by using Bisection\n";
                exit(1);
            }
            else {
                while (error > convergence_cri) {
                    fUpper=math::solvePolynomialsEq::subValToPolynomial(power,coefs,upper);
                    fLower=math::solvePolynomialsEq::subValToPolynomial(power,coefs,lower);
                    output = 0.5*(upper+lower);
                    fOutput=math::solvePolynomialsEq::subValToPolynomial(power,coefs,output);
                    if (fUpper*fOutput<0)
                    {
                        lower=output;
                    }
                    else if (fLower*fOutput<0) {
                        upper=output;
                    }
                    error = fabs((output - outputOld)/output);
                    outputOld=output;
                }
            }
            return output;
        }

        double subValToPolynomial(std::vector<double> &power, std::vector<double> &coefs, double Value)
        {
            int polySize(static_cast<int>(power.size()));
            double fValue(0.0);
            for (int polyOrder = 0; polyOrder < polySize; polyOrder++) {
                fValue+=pow(Value,power[polyOrder])*coefs[polyOrder];
            }
            return fValue;
        }

        std::tuple<bool, double, double> polynominal2d(double A, double B, double C)
        {
            double delta(B*B-4*A*C);
            if (delta>0)
            {
                return std::make_tuple(
                            true,
                            (-B+pow(delta,0.5))/(2*A),
                            (-B-pow(delta,0.5))/(2*A));
            }
            else if (delta==0.0)
            {
                return std::make_tuple(
                            true,
                            -B/(2*A),
                            -B/(2*A));
            }
            else
            {
                return std::make_tuple(false,0.0,0.0);
            }
        }
    }

    namespace BR2Fncs {
    double pointAuxValue_vol(int element, double a, double b, int valType, int dir)
    {
        double out(0.0);
        std::vector<double> Value(mathVar::orderElem + 1, 0.0);

        Value = auxUlti::getElementAuxValuesOfOrder_BR2_vol(element, valType, dir);

        math::basisFc(a, b, auxUlti::checkType(element));
        for (int order = 0; order <= mathVar::orderElem; order++)
        {
            out += Value[order] * mathVar::B[order];
        }
        return out;
    }

    double pointAuxValue_sur(int edge, int element, double a, double b, int valType, int dir)
    {
        double out(0.0);
        std::vector<double> Value(mathVar::orderElem + 1, 0.0);

        Value = auxUlti::getElementAuxValuesOfOrder_BR2_sur(edge, element, valType, dir);

        math::basisFc(a, b, auxUlti::checkType(element));
        for (int order = 0; order <= mathVar::orderElem; order++)
        {
            out += Value[order] * mathVar::B[order];
        }
        return out;
    }
    }

    double plusSideSurfaceValue(int edge, int element, int nG, int valType, int valKind)
    {
        double a(0.0), b(0.0);
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        return math::pointValue(element, a, b, valType, valKind);
    }

    double plusSideSurfaceDerivativeValue(int edge, int element, int nG, int valType, int dir)
    {
        /* ValType la id cua bien can tinh:
         * 1: (mu)dRho
         * 2: (mu)dRhou
         * 3: (mu)dRhov
         * 4: (mu)dRhoE
         * 5: dRho --> day la bien phu khi mass diffusion on: Sm = div(rho)
        */

        double valPlus(0.0), a(0.0), b(0.0);
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        if (systemVar::auxVariables==1)
        {
            valPlus = math::pointAuxValue(element, a, b, valType, dir);
        }
        else if (systemVar::auxVariables==2)
        {

        }

        return valPlus;
    }
}//end of namespace math
