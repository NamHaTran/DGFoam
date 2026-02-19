/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgRefCell.H"
#include "basisFunctions.H"
#include "error.H"
#include "vector.H"
#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dgRefCell::dgRefCell(const label pOrder, const dgCellType type)
:
    pOrder_(pOrder),
    type_(type)
{
    generateCellGaussPointsAndWeights();
    computeBasisAndDerivatives();
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::dgRefCell::getGaussLegendre1D
(
    const label pOrder,
    List<scalar>& eta,
    List<scalar>& w
)
{
    // --------------------------------------------------
    // Number of Gauss points
    // nGauss = pOrder + 1
    // --------------------------------------------------

    const label nGauss = pOrder + 1;

    eta.setSize(nGauss);
    w.setSize(nGauss);

    switch (nGauss)
    {
        case 1:
            eta[0] = 0.0;
            w  [0] = 2.0;
            break;

        case 2:
            eta = { -0.5773502691896257,
                     0.5773502691896257 };
            w   = {  1.0, 1.0 };
            break;

        case 3:
            eta = { -0.7745966692414834,
                     0.0,
                     0.7745966692414834 };
            w   = {  0.5555555555555556,
                     0.8888888888888888,
                     0.5555555555555556 };
            break;

        case 4:
            eta = { -0.8611363115940526,
                    -0.3399810435848563,
                     0.3399810435848563,
                     0.8611363115940526 };
            w   = {  0.3478548451374539,
                     0.6521451548625461,
                     0.6521451548625461,
                     0.3478548451374539 };
            break;

        case 5:
            eta = { -0.9061798459386640,
                    -0.5384693101056831,
                     0.0,
                     0.5384693101056831,
                     0.9061798459386640 };
            w   = {  0.2369268850561891,
                     0.4786286704993665,
                     0.5688888888888889,
                     0.4786286704993665,
                     0.2369268850561891 };
            break;

        case 6:
            eta = { -0.9324695142031521,
                    -0.6612093864662645,
                    -0.2386191860831969,
                     0.2386191860831969,
                     0.6612093864662645,
                     0.9324695142031521 };
            w   = {  0.1713244923791704,
                     0.3607615730481386,
                     0.4679139345726910,
                     0.4679139345726910,
                     0.3607615730481386,
                     0.1713244923791704 };
            break;

        case 7:
            eta = { -0.9491079123427585,
                    -0.7415311855993945,
                    -0.4058451513773972,
                     0.0,
                     0.4058451513773972,
                     0.7415311855993945,
                     0.9491079123427585 };
            w   = {  0.1294849661688697,
                     0.2797053914892766,
                     0.3818300505051189,
                     0.4179591836734694,
                     0.3818300505051189,
                     0.2797053914892766,
                     0.1294849661688697 };
            break;

        case 8:
            eta = { -0.9602898564975363,
                    -0.7966664774136267,
                    -0.5255324099163290,
                    -0.1834346424956498,
                     0.1834346424956498,
                     0.5255324099163290,
                     0.7966664774136267,
                     0.9602898564975363 };
            w   = {  0.1012285362903763,
                     0.2223810344533745,
                     0.3137066458778873,
                     0.3626837833783620,
                     0.3626837833783620,
                     0.3137066458778873,
                     0.2223810344533745,
                     0.1012285362903763 };
            break;

        case 9:
            eta = { -0.9681602395076261,
                    -0.8360311073266358,
                    -0.6133714327005904,
                    -0.3242534234038089,
                     0.0,
                     0.3242534234038089,
                     0.6133714327005904,
                     0.8360311073266358,
                     0.9681602395076261 };
            w   = {  0.0812743883615744,
                     0.1806481606948574,
                     0.2606106964029354,
                     0.3123470770400029,
                     0.3302393550012598,
                     0.3123470770400029,
                     0.2606106964029354,
                     0.1806481606948574,
                     0.0812743883615744 };
            break;

        case 10:
            eta = { -0.9739065285171717,
                    -0.8650633666889845,
                    -0.6794095682990244,
                    -0.4333953941292472,
                    -0.1488743389816312,
                     0.1488743389816312,
                     0.4333953941292472,
                     0.6794095682990244,
                     0.8650633666889845,
                     0.9739065285171717 };
            w   = {  0.0666713443086881,
                     0.1494513491505806,
                     0.2190863625159820,
                     0.2692667193099963,
                     0.2955242247147529,
                     0.2955242247147529,
                     0.2692667193099963,
                     0.2190863625159820,
                     0.1494513491505806,
                     0.0666713443086881 };
            break;

        default:
            FatalErrorInFunction
                << "Gauss-Legendre not implemented for nGauss = "
                << nGauss << nl
                << "Derived from pOrder = " << pOrder
                << abort(FatalError);
    }
}

/*
void Foam::dgRefCell::getGaussJacobi1D
(
    const label pOrder,
    List<scalar>& eta,
    List<scalar>& w
)
{
    // --------------------------------------------------
    // Fallback conditions
    // --------------------------------------------------
    if (pOrder == 0 || pOrder > 4)
    {
        getGaussLegendre1D(pOrder, eta, w);
        return;
    }

    const label nGauss = pOrder + 1;
    eta.setSize(nGauss);
    w.setSize(nGauss);

    switch (nGauss)
    {
        case 2:
            eta = { -0.6546536707079771,
                    0.1546536707079771 };
            w  = {  0.7111111111111111,
                    0.6222222222222222 };
            break;

        case 3:
            eta = { -0.7660444431189780,
                   -0.1736481776669303,
                    0.5 };
            w  = {  0.4444444444444444,
                    0.7777777777777778,
                    0.1111111111111111 };
            break;

        case 4:
            eta = { -0.8302238962785669,
                   -0.4688487934707142,
                    0.0558466708154264,
                    0.7168821974518675 };
            w  = {  0.3038329932404390,
                    0.5909114117464944,
                    0.3719183942939050,
                    0.0666708673851616 };
            break;

        default:
            FatalErrorInFunction
                << "Invalid nGauss for Gauss-Jacobi"
                << abort(FatalError);
    }
}
*/


void Foam::dgRefCell::getGaussLobatto1D
(
    const label pOrder,
    List<scalar>& eta,
    List<scalar>& w
)
{
    // --------------------------------------------------
    // Number of GLL points = pOrder + 2
    // --------------------------------------------------
    const label nGauss = pOrder + 2;

    // --------------------------------------------------
    // Supported range: 2 <= nGauss <= 10
    // If outside → fallback to Gauss-Legendre
    // --------------------------------------------------
    if (nGauss < 2 || nGauss > 10)
    {
        getGaussLegendre1D(pOrder, eta, w);
        return;
    }

    eta.setSize(nGauss);
    w.setSize(nGauss);

    switch (nGauss)
    {
        case 2:
            eta = { -1.0, 1.0 };
            w   = {  1.0, 1.0 };
            break;

        case 3:
            eta = { -1.0, 0.0, 1.0 };
            w   = {  0.3333333333333333,
                     1.3333333333333333,
                     0.3333333333333333 };
            break;

        case 4:
            eta = { -1.0,
                   -0.4472135954999579,
                    0.4472135954999579,
                    1.0 };
            w   = {  0.1666666666666667,
                     0.8333333333333333,
                     0.8333333333333333,
                     0.1666666666666667 };
            break;

        case 5:
            eta = { -1.0,
                   -0.6546536707079771,
                    0.0,
                    0.6546536707079771,
                    1.0 };
            w   = {  0.1,
                     0.5444444444444444,
                     0.7111111111111111,
                     0.5444444444444444,
                     0.1 };
            break;

        case 6:
            eta = { -1.0,
                   -0.7650553239294647,
                   -0.2852315164806451,
                    0.2852315164806451,
                    0.7650553239294647,
                    1.0 };
            w   = {  0.0666666666666667,
                     0.3784749562978470,
                     0.5548583770354864,
                     0.5548583770354864,
                     0.3784749562978470,
                     0.0666666666666667 };
            break;

        case 7:
            eta = { -1.0,
                   -0.8302238962785670,
                   -0.4688487934707142,
                    0.0,
                    0.4688487934707142,
                    0.8302238962785670,
                    1.0 };
            w   = {  0.0476190476190476,
                     0.2768260473615659,
                     0.4317453812098626,
                     0.4876190476190476,
                     0.4317453812098626,
                     0.2768260473615659,
                     0.0476190476190476 };
            break;

        case 8:
            eta = { -1.0,
                   -0.8717401485096066,
                   -0.5917001814331423,
                   -0.2092992179024789,
                    0.2092992179024789,
                    0.5917001814331423,
                    0.8717401485096066,
                    1.0 };
            w   = {  0.0357142857142857,
                     0.2107042271435061,
                     0.3411226924835044,
                     0.4124587946587039,
                     0.4124587946587039,
                     0.3411226924835044,
                     0.2107042271435061,
                     0.0357142857142857 };
            break;

        case 9:
            eta = { -1.0,
                   -0.8997579954114602,
                   -0.6771862795107377,
                   -0.3631174638261782,
                    0.0,
                    0.3631174638261782,
                    0.6771862795107377,
                    0.8997579954114602,
                    1.0 };
            w   = {  0.0277777777777778,
                     0.1654953615608055,
                     0.2745387125001617,
                     0.3464285109730463,
                     0.3715192743764172,
                     0.3464285109730463,
                     0.2745387125001617,
                     0.1654953615608055,
                     0.0277777777777778 };
            break;

        case 10:
            eta = { -1.0,
                   -0.9195339081664589,
                   -0.7387738651055051,
                   -0.4779249498104445,
                   -0.1652789576663870,
                    0.1652789576663870,
                    0.4779249498104445,
                    0.7387738651055051,
                    0.9195339081664589,
                    1.0 };
            w   = {  0.0222222222222222,
                     0.1333059908510701,
                     0.2248893420631265,
                     0.2920426836796838,
                     0.3275397611838975,
                     0.3275397611838975,
                     0.2920426836796838,
                     0.2248893420631265,
                     0.1333059908510701,
                     0.0222222222222222 };
            break;

        default:
            FatalErrorInFunction
                << "Invalid nGauss for Gauss-Lobatto"
                << abort(FatalError);
    }
}

void Foam::dgRefCell::getGaussRadau1D
(
    const label pOrder,
    const bool includeLeft,   // true  -> include eta = -1  (Left Radau)
                              // false -> include eta = +1  (Right Radau)
    List<scalar>& eta,
    List<scalar>& w
)
{
    // --------------------------------------------------
    // Number of quadrature points
    // Radau: n = pOrder + 1
    // --------------------------------------------------
    const label nGauss = pOrder + 2;

    // --------------------------------------------------
    // Supported range: 2 <= nGauss <= 10
    // --------------------------------------------------
    if (nGauss < 2 || nGauss > 10)
    {
        getGaussLegendre1D(pOrder, eta, w);
        return;
    }

    eta.setSize(nGauss);
    w.setSize(nGauss);

    // --------------------------------------------------
    // Left Radau (include eta = -1)
    // --------------------------------------------------
    if (includeLeft)
    {
        switch (nGauss)
        {
            case 2:
                eta = { -1.0,  0.3333333333333333 };
                w   = {  0.5,  1.5 };
                break;

            case 3:
                eta = { -1.0,
                        -0.2898979485566356,
                         0.6898979485566356 };
                w   = {  0.2222222222222222,
                         1.0249716523768432,
                         0.7528061254009346 };
                break;

            case 4:
                eta = { -1.0,
                        -0.5753189235216941,
                         0.1810662711185306,
                         0.8228240809745921 };
                w   = {  0.1250000000000000,
                         0.6576886399601195,
                         0.7763869376863437,
                         0.4409244223535368 };
                break;

            case 5:
                eta = { -1.0,
                        -0.7204802713124389,
                        -0.1671808647378336,
                         0.4463139727237523,
                         0.8857916077709646 };
                w   = {  0.0800000000000000,
                         0.4462078021671415,
                         0.6236530459514825,
                         0.5627120302989241,
                         0.2874271215824519 };
                break;

            case 6:
                eta = { -1.0,
                        -0.8029298284023471,
                        -0.3909285467072722,
                         0.1240503795052277,
                         0.6039731642527837,
                         0.9203802858970625 };
                w   = {  0.0555555555555556,
                         0.3196407532205109,
                         0.4853871884689699,
                         0.5209267831895740,
                         0.4169013343119077,
                         0.2015883852534816 };
                break;

            case 7:
                eta = { -1.0,
                        -0.8538913426394822,
                        -0.5384677240601090,
                        -0.1173430375435740,
                         0.3260306194376914,
                         0.7038428006630314,
                         0.9413671456804302 };
                w   = {  0.0408163265306122,
                         0.2392274892253124,
                         0.3809498736442311,
                         0.4471098290145665,
                         0.4247037790059556,
                         0.3182042314679501,
                         0.1489884711119649 };
                break;

            case 8:
                eta = { -1.0,
                        -0.8874748789261557,
                        -0.6395186165269650,
                        -0.2947505657736607,
                         0.0943072526611108,
                         0.4684203544308211,
                         0.7706418936781916,
                         0.9550412271225750 };
                w   = {  0.0312500000000000,
                         0.1853581548029793,
                         0.3041306206467851,
                         0.3765175453891186,
                         0.3915721674524936,
                         0.3470147956345014,
                         0.2496479013298649,
                         0.1145088147442560 };
                break;

            case 9:
                eta = { -1.0,
                        -0.9107320894200603,
                        -0.7112674859157086,
                        -0.4263504857111386,
                        -0.0903733696068530,
                         0.2561356708334554,
                         0.5713830412087385,
                         0.8173527842004121,
                         0.9644401697052731 };
                w   = {  0.0246913580246914,
                         0.1476540190469377,
                         0.2471893782045939,
                         0.3168437756704379,
                         0.3482729979163132,
                         0.3376939669759299,
                         0.2863866963572311,
                         0.2005532989972816,
                         0.0907145088065986 };
                break;

            case 10:
                eta = { -1.0,
                        -0.9274843742335811,
                        -0.7638420424200026,
                        -0.5256460303700790,
                        -0.2362344693905880,
                         0.0760591978379781,
                         0.3806648401447244,
                         0.6477666876740094,
                         0.8512252205816079,
                         0.9711751807022469 };
                w   = {  0.0200000000000000,
                         0.1202966705574816,
                         0.2042701318799835,
                         0.2681948378419016,
                         0.3058592877244226,
                         0.3135824572269384,
                         0.2906101648323127,
                         0.2391934317143797,
                         0.1643760127369214,
                         0.0736170048456588 };
                break;

            default:
                FatalErrorInFunction
                    << "Invalid nGauss for Gauss-Radau rule"
                    << abort(FatalError);
        }
    }
    // --------------------------------------------------
    // Right Radau (mirror of Left Radau)
    // --------------------------------------------------
    else
    {
        List<scalar> etaLeft, wLeft;

        // Recursive call to generate left version
        getGaussRadau1D(pOrder, true, etaLeft, wLeft);

        // Mirror: x_i -> -x_{n-1-i}
        for (label i = 0; i < nGauss; ++i)
        {
            eta[i] = -etaLeft[nGauss - 1 - i];
            w[i]   =  wLeft[nGauss - 1 - i];
        }
    }
}

void Foam::dgRefCell::generateCellGaussPointsAndWeights()
{
    // --------------------------------------------------
    // 1D Gauss points and weights for each direction
    // --------------------------------------------------
    List<scalar> eta1, w1;
    List<scalar> eta2, w2;
    List<scalar> eta3, w3;

    // --------------------------------------------------
    // Select quadrature rules based on cell type
    // --------------------------------------------------
    switch (type_)
    {
        // ==================================================
        // Hexahedral cell
        // Tensor-product Gauss–Legendre in all directions
        // ==================================================
        case dgCellType::HEX:
        {
            getGaussLobatto1D(pOrder_, eta1, w1);
            getGaussLobatto1D(pOrder_, eta2, w2);
            getGaussLobatto1D(pOrder_, eta3, w3);
            break;
        }

        // ==================================================
        // Tetrahedral cell
        // Collapsed coordinates:
        //   - eta1: Gauss–Jacobi (alpha=2, beta=0)
        //   - eta2, eta3: Gauss–Legendre
        // pOrder = 0 handled internally by fallback
        // ==================================================
        case dgCellType::TET:
        {
            getGaussLobatto1D(pOrder_, eta1, w1);
            getGaussRadau1D(pOrder_, true, eta2, w2);
            getGaussRadau1D(pOrder_, true, eta3, w3);
            break;
        }

        // ==================================================
        // Prism cell
        // Collapsed direction in eta3:
        //   - eta1, eta2: Gauss–Lobatto
        //   - eta3: Gauss–Radau
        // pOrder = 0 handled internally by fallback
        // ==================================================
        case dgCellType::PRISM:
        {
            getGaussLobatto1D(pOrder_, eta1, w1);
            getGaussLobatto1D(pOrder_, eta2, w2);
            getGaussRadau1D  (pOrder_, true, eta3, w3);
            break;
        }

        // ==================================================
        // Pyramid cell
        // Same quadrature strategy as prism
        // ==================================================
        case dgCellType::PYRAMID:
        {
            getGaussLobatto1D(pOrder_, eta1, w1);
            getGaussLobatto1D(pOrder_, eta2, w2);
            getGaussRadau1D  (pOrder_, true, eta3, w3);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unsupported dgCellType"
                << abort(FatalError);
        }
    }

    // --------------------------------------------------
    // Total number of Gauss points (tensor-product)
    // --------------------------------------------------
    nGauss_ = eta1.size()*eta2.size()*eta3.size();

    gaussP_.setSize(nGauss_);
    wGauss_.setSize(nGauss_);

    // --------------------------------------------------
    // Build tensor-product Gauss points and weights
    // Reference coordinates: (eta1, eta2, eta3)
    // --------------------------------------------------
    label idx = 0;

    for (label k = 0; k < eta3.size(); ++k)
    {
        for (label j = 0; j < eta2.size(); ++j)
        {
            for (label i = 0; i < eta1.size(); ++i)
            {
                gaussP_[idx] =
                    vector
                    (
                        eta1[i],
                        eta2[j],
                        eta3[k]
                    );

                wGauss_[idx] =
                    w1[i]*w2[j]*w3[k];

                ++idx;
            }
        }
    }
}


void dgRefCell::computeBasisAndDerivatives()
{
    const label nBasis = Foam::getNumBasis(pOrder_, type_);
    nDof_ = nBasis;

    basis_.setSize(nGauss_);
    dBasis_dEta1_.setSize(nGauss_);
    dBasis_dEta2_.setSize(nGauss_);
    dBasis_dEta3_.setSize(nGauss_);

    for (label gp = 0; gp < nGauss_; ++gp)
    {
        switch (type_)
        {
            case dgCellType::HEX:
                Foam::computeHexBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PRISM:
                Foam::computePrismBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::TET:
                Foam::computeTetBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            case dgCellType::PYRAMID:
                Foam::computePyramidBasisAndDerivatives(
                    gaussP_[gp], pOrder_,
                    basis_[gp],
                    dBasis_dEta1_[gp],
                    dBasis_dEta2_[gp],
                    dBasis_dEta3_[gp]
                );
                break;

            default:
                FatalErrorInFunction
                    << "Unsupported cell type: " << type_ << nl
                    << abort(FatalError);
        }
    }
}

// ************************************************************************* //
